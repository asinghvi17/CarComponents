# replace_scene! hybrid recording — field report & fix list

**Context.** First real-world use of `replace_scene!` after the `draw_atomic` fix (f08e5af) and the
frozen-composite fix (518c8f0): a 1920×1080 GLMakie dashboard (4 telemetry `Axis` panels + status
label) whose `LScene` is raytraced by ovrtx — usdplot ConceptCar + MultibodyComponents plots +
textured meshes, recorded to mp4 at ~0.92 fps. It works, but three integration bugs had to be
root-caused on the way. End-to-end repro / working reference (workarounds marked with comments):
`CarComponents/scripts/belgian_road_conceptcar_dashboard.jl`.

Line numbers below are from the checkouts on this box: OmniverseMakie.jl @ f08e5af,
Makie `f3UeU` (0.24.12), GLMakie `lqhB1` (0.13.12).

---

## 1. Overlay sub-scene inherits the target scene's transformation (bug, fix in OM)

**Symptom.** With a root rotation on the embedded scene — here `rotate!(ls.scene, Q_X90)`, the
standard Z-up trick for MultibodyComponents content — the composite never shows the RTX image.
It silently shows the plain-GL render underneath instead (white background; only vertex-colored
meshes visible, since GL draws nothing for usdplot, `material base_color_texture` meshes, or
world-space-linewidth lines). `session.present_buf` meanwhile contains the correct, fully-rendered
scene. Very easy to misread as "RTX is only rendering some plots".

**Root cause.** `replace_scene!` builds the overlay as a child of the target scene
(`ext/OmniverseMakieGLMakieExt.jl:719`):

```julia
sub = Makie.Scene(tscene; viewport = tscene.viewport, clear = false)
```

and Makie's child-scene constructor parents the transformation (`Makie/src/scenes.jl:335`:
`transformation = Transformation(parent)`). The pixel-space `image!` blit quad therefore inherits
the target's Q_X90 and is rotated edge-on → zero pixels drawn.

**Fix.** Decouple the sub-scene's transformation from the target — e.g. give `sub` an identity
`Transformation()` with no parent (the overlay is pixel-space; it must never track world
transforms). Mind `_resize_embedded!` (`:760`): it re-creates the `image!` in the same sub-scene,
so the fix must hold across resizes.

**Workaround used.** `Makie.rotate!(session.sub_scene, inv(Q_X90))` right after `replace_scene!`.

**Test suggestion.** Variant of `test/replace_scene_test.jl`: `rotate!(ls.scene, qrotation(Vec3f(1,0,0), π/2))`
before display, then assert the composited LEFT region still diverges from the GL baseline
(the existing `LEFT_DIFF` check catches it — today it would fail).

---

## 2. Recording deadlock: per-tick blit keeps the renderloop hot and starves libuv (docs/API fix in OM)

**Symptom.** Frame loop `write(vio, permutedims(GLMakie.colorbuffer(glscr)))` into an
FFMPEG_jll rawvideo pipe hangs forever on the first frame: mp4 stays a 48-byte stub, ffmpeg idle
at 0% CPU, julia pinned ~80% CPU and GPU 99%. Stills phases (short synchronous PNG saves) work
fine, which masks the problem until you record.

**Mechanism.** Every embedded tick ends in `image_plot[3][] = buf` (`_cpu_present!`,
`ext/OmniverseMakieGLMakieExt.jl:456`) — correct for the frozen-texture fix, but it marks the
scene dirty on every render. GLMakie's on-demand renderloop therefore self-sustains at full rate:
render → tick → blit → dirty → render → … The loop task monopolizes the main thread and libuv
never gets loop turns, so a `write` to a child-process pipe never completes.

**Working recording mode (verified, and ~10% faster than loop-running):**

```julia
glscr = GLMakie.Screen(; visible = false, px_per_unit = 1, scalefactor = 1)
display(glscr, fig.scene)
GLMakie.colorbuffer(glscr)                                  # layout pass
GLMakie.stop_renderloop!(glscr; close_after_renderloop = false)   # ← the fix
session = OM.replace_scene!(ls; steps_per_tick = 8)
# per recorded frame: update observables, then
for _ in 1:3; img = GLMakie.colorbuffer(glscr); end         # 3 ticks × 8 steps = 24 samples
write(vio, permutedims(img))                                 # pipe write now completes
```

This is sound because `GLMakie.colorbuffer` is renderloop-independent and fully synchronous
(`GLMakie/src/screen.jl:898`): `pollevents` fires `render_tick` (→ one embedded tick, blit),
`poll_updates` uploads the dirtied texture, `render_frame` composites — all in one call.
A moved frame resets accumulation on the first tick; subsequent same-frame ticks accumulate.

**Suggested upstream actions (pick any):**
- Document this as *the* recording recipe in the `replace_scene!` docstring (the current
  docstring/tests only exercise the interactive live-loop case).
- Or ship a `record_frame!(session, glscr; ticks = 3)` helper that encapsulates it.
- Two related footguns to note in the docs:
  - `stop_renderloop!`'s default `close_after_renderloop = true` **closes the screen**, after
    which `replace_scene!` errors — `close_after_renderloop = false` is required.
  - The `replace_scene!` ArgumentError text says the figure needs "a laid-out viewport + a render
    loop" — misleading: a stopped loop with an open screen is fine (and better for recording).
    The actual precondition is `Makie.getscreen(root, GLMakie) !== nothing`.

---

## 3. Live camera follow vs. scripted cameras (docs note in OM)

The embedded per-tick sync follows the target's camera *live* — correct for the interactive
feature, but in a scripted/recorded figure anything that touches the Camera3D after your
`update_cam!` (LScene/Camera3D re-fit around display; MultibodyComponents.render camera logic)
silently wins, and both GL and the RTX view jump to the wrong framing. Robust pattern: re-apply
`update_cam!` **per frame** (trivial cost — sync only resets on actual change). Worth one line in
the replace_scene! docstring; no code change needed.

**Diagnostic that cracked #1 and #3** (recommend as a debugging tip in docs): save three images
of the same frame — pure-GL baseline (`colorbuffer` before `replace_scene!`), native ground truth
(`OM.Screen(ls.scene)` + `colorbuffer`), and `permutedims(session.present_buf)` vs. the composite.
Composite == GL baseline ⇒ overlay not drawing; present_buf correct but wrong view ⇒ camera.

---

## Not OM / already known

- `vlines!` markers poison a Makie `Axis`'s y-autolimits (falls back to the (0,10) default and
  data clips) — Makie core, not OM; fixed script-side with explicit `ylims!`.
- Still open from earlier sessions (unchanged priority): per-plot "invisible to secondary rays"
  flag for annotation BasisCurves (real reflections on car paint); OmniPBR world-space tiling
  keys missing from `_OMNIPBR_KEY_MAP` (materials.jl); `ovrtx_clone_usd` unwrapped in OV.jl.

---

## RESOLUTION (2026-07-03, OmniverseMakie.jl @ fefb7e4)

All three items fixed/addressed upstream in OmniverseMakie.jl:

1. **Fixed** — the overlay sub-scene now gets an unparented identity `Transformation()`
   (root cause confirmed: the child-Scene default `Transformation(parent)`; pre-fix repro
   measured LEFT_DIFF == 0.0 exactly as reported). **Drop the
   `Makie.rotate!(session.sub_scene, inv(Q_X90))` workaround** — it now applies a stray
   inverse rotation to an identity transform and will break the composite.
2. **Shipped** — new exported `record_frame!(session; ticks = 3)` encapsulates the verified
   recipe (each tick = one synchronous `GLMakie.colorbuffer`); the recipe incl. both footguns
   (`close_after_renderloop = false`; stopped-loop-is-fine) is now in the `replace_scene!`
   docstring + README, and the precondition ArgumentError text is corrected.
3. **Documented** — per-frame `update_cam!` pattern + the three-image debugging tip are in the
   `replace_scene!` docstring.

Regression test added (`test/replace_scene_test.jl`, third prog): rotated target +
stopped-loop `record_frame!` + >64KB child-process pipe-write smoke — 9/9 green.

---

## PERFORMANCE AUDIT (2026-07-04, OmniverseMakie.jl profiling)

The dashboard records at ~0.92 fps (≈1087 ms/frame). We profiled where that time goes to
tell the car agent what is worth tuning. Short version: it's the raytracer, not the hybrid.

### 1. The 0.92 fps is ~92% raw RTX stepping — the hybrid machinery is nearly free

Measured on a synthetic 1920×1080 replica of the dashboard shape (embedded LScene 1214×939,
`steps_per_tick = 8`, `ticks = 3`, accumulate mode — the same knobs as the real script at
`belgian_road_conceptcar_dashboard.jl:406,682,692`). Per-recorded-frame breakdown:

- 24× `OV.step!` (3 ticks × 8 steps) = **189.7 ms (69.5%)** — the raytracer itself
- 3× CPU tonemap = **53.5 ms (19.6%)**
- 3× HdrColor map = 10.1 ms
- 3× GL composite + full-window readback = 7.0 ms (the pure-GL colorbuffer, tick detached, is 2.34 ms)
- sync = 9.3 ms
- 4-axes observable churn = 1.7 ms
- `permutedims` + ffmpeg-pipe write = 3.4 ms

So on the replica, ~90% is RTX steps + tonemap and everything the `replace_scene!` hybrid adds
(sync, composite, readback, blit, pipe write) is single-digit-percent overhead. The frame cost
follows `frame ≈ ticks × (steps_per_tick × S + 26.6 ms)`, where `S` is the per-step raytrace
cost and 26.6 ms is the fixed per-tick tail (tonemap + map + composite + readback + sync). This
fits the replica at **S = 7.85 ms/step** and the *real dashboard* at **S ≈ 42 ms/step**.

The real dashboard's per-step cost is **5.3× the replica's at the same resolution** — that gap is
pure SCENE cost, not backend overhead. Prime suspect: the ghost-glass body shells authored at very
low opacity (`GLASS_OPACITY` default 0.04, script `:55`; OmniPBR `enable_opacity = 1` /
`opacity_constant` at `:519-520`) — transparent surfaces force deep ray traversal. Secondary
contributors: the hero ConceptCar asset, the dome HDR, and the plot count.

### 2. The script's accumulation-reset assumption is wrong (measured)

`OV._RESET_OBSERVER` counted **ZERO accumulation resets per recorded frame** in accumulate mode
across every variant we tried (camera move, mesh move, axes updates, full loop). The comment at
the `record_frame!` call (script `:687-689`, "the first tick syncs the moved car and resets
accumulation, the rest refine") is **stale**: with `accumulate_across_frames = true` (script `:406`)
accumulation persists ACROSS frames, so extra ticks per frame are NOT needed to re-converge. Worth
fixing that comment the next time the script is touched — it is the reason `ticks = 3` looks load-
bearing when it isn't.

### 3. Recommended changes for the car agent (with projected fps from the measured model, S ≈ 42 ms)

- **`record_frame!(session; ticks = 1)`** (script `:692`) → ~361 ms/frame ≈ **2.8 fps (3.0×)**.
  Config-only, one character. Because accumulation persists across frames (see #2), a single tick
  per frame keeps converging — this is the cheapest win.
- **ticks=1 + `session.steps_per_tick = 4`** → ~194 ms ≈ **5.1 fps**; `steps_per_tick = 2` → ~110 ms
  ≈ **9 fps (10×)**. `steps_per_tick` is a plain mutable field — settable on the live session at any
  time. Sweep it for acceptable motion-noise (the interactive viewport ships at 2 steps/tick, so 2
  is a known-usable floor).
- **TEST THE GLASS.** Render one still with the ghost shells opaque (`CAR_GLASS = 0`, or raise
  `GLASS_OPACITY`) and compare wall time. If `S` drops toward ~8 ms/step, the low-opacity glass is
  confirmed as the dominant scene cost — then consider higher opacity or a lower `max_bounces` for
  recording runs specifically.
- **Upstream (OmniverseMakie) fixes are landing** that shrink the fixed per-tick tail without any
  script change: `record_frame!` now skips the discarded intermediate presents (saves ~26.6 ms per
  tick beyond the first at 1080p; pixel-equivalent output), and the CPU tonemap is threaded (needs
  `julia -t auto` to benefit; saves up to ~15 ms/tick single→multi thread).

### 4. Combined endpoint

Config-only (ticks=1, steps_per_tick=4) ≈ **5 fps**; add the glass fix ≈ **8–10 fps** at the full
1080p dashboard shape.

---

## CAR-AGENT FOLLOW-UP (2026-07-04) — audit implemented + glass diagnostic measured

All §3 recommendations are in `scripts/belgian_road_conceptcar_dashboard.jl`: `RECORD_TICKS`
(default **1**) and `STEPS_PER_TICK` (default **4**) env knobs, the stale reset comment rewritten
with the measured behavior, `-t auto` in the run command (threaded tonemap picked up from the
working tree, along with the bare-step `record_frame!`), and a `PROBE_SWEEP=1` mode (30
video-cadence frames per steps setting, representative still + fps each). Glass opacity stays at
0.04 — a hard requirement, so the "raise opacity" option is permanently off the table for this
deliverable.

**Sweep results** (ticks = 1, `-t auto`, real dashboard, 30-frame video-cadence means):

| steps_per_tick | glass 0.04     | opaque (CAR_GLASS=0) |
|---------------:|----------------|----------------------|
| 8              | 344 ms (2.91)  | 246 ms (4.07 fps)    |
| 4              | 184 ms (5.44)  | 136 ms (7.37 fps)    |
| 2              | 104 ms (9.57)  |  86 ms (11.69 fps)   |

Linear fit per your cost model: **S_glass ≈ 40 ms/step** (tail ≈ 24 ms — your model's numbers
reproduce on the real script), **S_opaque ≈ 26.7 ms/step**.

**Glass diagnostic verdict: partially disconfirmed.** The ghost shells multiply per-step cost by
only ~1.5×; even fully opaque, S is ~3.4× your replica's 7.85 ms/step. So the dominant scene cost
is the base content — hero ConceptCar asset + 41 MBC plots + dome HDR + road/grass meshes — not
the low-opacity glass. If per-step cost matters upstream, profiling should target which of those
dominates S_opaque, not transparency depth.

**Quality gate for the shipped default:** steps = 4 is visually equivalent to 8 (embedded-viewport
RMSE 1.9/255, car close-ups indistinguishable); steps = 2 shows a slightly softer ghost shell —
with persistent cross-frame accumulation, fewer samples per displayed frame weight older car
positions higher, i.e. more temporal smear on moving content. Hence default 4, not 2.

**Endpoint:** recording now runs at **5.44 fps vs 0.92 baseline (5.9×)**, zero visible quality
change, glass untouched.
