# lib/OpenCRG/src/transform.jl

"""
    integrate_reference_line(refline, phi) -> (x, y)

Integrate heading `phi` (`phi[i]` = heading of the segment ARRIVING at node
`i`; `phi[1]` is unused as an arrival heading — it's `REFERENCE_LINE_START_PHI`,
already substituted in by `assemble_channels`) into world-frame positions,
replicating `calcRefLine` in the reference implementation's `crgLoader.c`:

- No end position given (`refline.end_x === nothing` AND `refline.end_y ===
  nothing`): simple forward Euler, `x[i+1] = x[i] + du*cos(phi[i+1])`.
- Both `end_x`/`end_y` given: integrate backward from the true end point
  using the same arrival-phi convention, then blend the forward-continuation
  value with the backward-derived value linearly across the whole line
  (fraction 0 at the start, 1 at the end) — redistributing integration error
  instead of leaving it all as a discontinuity at the very end.
- Exactly one of `end_x`/`end_y` given: an error, not a silent fallback to
  Case A -- you can't sensibly blend only one of two coordinates, and a real
  file with only one of the pair set is far more likely to indicate a
  parsing bug or corruption than genuine spec-compliant data (the ASAM
  example files always declare both or neither).

`nu == 1` degenerates safely in both cases: the loop ranges become empty
(`1:0`), so only the start point is ever produced, and `end_x`/`end_y` are
never read.
"""
function integrate_reference_line(refline::ReferenceLineParams, phi::Vector{Float64})
    nu = length(phi)
    du = refline.increment
    x, y = Vector{Float64}(undef, nu), Vector{Float64}(undef, nu)
    x[1], y[1] = refline.start_x, refline.start_y
    has_end_x, has_end_y = refline.end_x !== nothing, refline.end_y !== nothing
    xor(has_end_x, has_end_y) &&
        error("REFERENCE_LINE_END_X/END_Y must both be set or both be absent, not just one")
    if !has_end_x && !has_end_y
        for i in 1:nu-1
            x[i+1] = x[i] + du * cos(phi[i+1])
            y[i+1] = y[i] + du * sin(phi[i+1])
        end
        return x, y
    end
    xb, yb = Vector{Float64}(undef, nu), Vector{Float64}(undef, nu)
    xb[nu], yb[nu] = refline.end_x, refline.end_y
    for i in nu:-1:2
        xb[i-1] = xb[i] - du * cos(phi[i])
        yb[i-1] = yb[i] - du * sin(phi[i])
    end
    for i in 1:nu-1
        fraction = i / (nu - 1)
        x[i+1] = (1 - fraction) * (x[i] + du * cos(phi[i+1])) + fraction * xb[i+1]
        y[i+1] = (1 - fraction) * (y[i] + du * sin(phi[i+1])) + fraction * yb[i+1]
    end
    return x, y
end

"""
    integrate_reference_z(refline, slope, nu) -> z_ref

1D analogue of `integrate_reference_line`, integrating longitudinal `slope`
(dz/du) into a reference-line elevation profile. If there's no `slope`
channel at all, `refline.start_slope` (default 0.0, from
`REFERENCE_LINE_START_S`) is used as a constant slope for every step. If
there's no slope channel AND `start_slope == 0.0`, the reference
implementation (`calcRefLineZ`) skips integration early — but this does NOT
mean elevation is zero everywhere. Cross-checked against the real compiled
C reference library (`crgEvaluv2z`/`crgEvalz.c`): when the ref-z channel is
invalid, it falls back to adding `channelRefZ.info.first` (i.e.
`REFERENCE_LINE_START_Z`) at every point, not zero. So a flat road at a
nonzero elevation (`start_z` set, no slope) must produce a constant
`start_z` vector, not a constant-zero one — `end_z` is irrelevant in this
branch (confirmed against the oracle: it's ignored even when set), since
there's no slope data to integrate toward it anyway.

`nu == 1` degenerates safely (only `z_ref[1] = refline.start_z` is ever
produced; `end_z` is unused), matching `integrate_reference_line`.

**End-anchored blend fraction is deliberately NOT `i/(nu-1)` by analogy with
`integrate_reference_line`'s (x,y) blend — found and fixed during Task 17's
cross-validation.** `integrate_reference_line`'s blend (`calcRefLine` in
`crgLoader.c`, ~line 2181-2190) and this function's blend
(`calcRefLineZ`, ~line 2260-2268) look structurally parallel but are NOT
numerically identical in the reference C implementation: `calcRefLine` uses
`fraction = (i+1)/(size-1)` and loops `i` from `0` to `size-2` inclusive
(so every index from 1 to `size-1`, i.e. including the anchor, gets
"re-blended" — a no-op there since `fraction=1` at the last step reproduces
the anchor exactly); `calcRefLineZ` instead uses `fraction = i/(size-1)`
(the loop variable itself, NOT `i+1`) and loops only up to `size-3`
inclusive (`i < size-2`), one iteration short of `calcRefLine`'s bound — so
the very last index (`size-1`) is never revisited by the blend formula at
all; it keeps the value written directly at the top of the function
(the true `end_z` anchor). Net effect: at a given target node, `calcRefLineZ`
uses a fraction ONE NODE "behind" what the symmetric-by-analogy formula
would give. This was invisible to every pre-Task-17 test (both
`integrate_reference_z`'s own unit tests and every cross-validation fixture)
because the only fixture that exercises this branch at all
(`synthetic_end_anchored.crg`, Task 11) declares just one long-section
column, which the real C library refuses to load at all (`crgLoader.c`:
"no or insufficient long section data available", requires >= 2) — so this
branch was cross-validated against pure hand arithmetic only, never against
the compiled oracle, until Task 17 added `synthetic_end_anchored_2col.crg`.
Cross-validating that new fixture's z values against `crgEvaluv2z` surfaced
a ~0.1-0.15m systematic mismatch at every interior node, which traced
directly to this fraction/loop-bound discrepancy (verified by hand against
the disassembled formula above, and confirmed bit-for-bit against the C
oracle's actual `channelRefZ.data[]` values once ported). Whether this
asymmetry in the C reference is an intentional design choice or an
unnoticed bug upstream is unknown and irrelevant here — cross-validating
against the compiled oracle, not against what "should" be symmetric, is
this whole package's stated design philosophy, so this function replicates
`calcRefLineZ` exactly rather than `calcRefLine`'s pattern.
"""
function integrate_reference_z(refline::ReferenceLineParams, slope::Union{Vector{Float64},Nothing}, nu::Int)
    if slope === nothing && refline.start_slope == 0.0
        return fill(refline.start_z, nu)
    end
    du = refline.increment
    slope_at(i) = slope === nothing ? refline.start_slope : slope[i]
    z_ref = Vector{Float64}(undef, nu)
    z_ref[1] = refline.start_z
    if refline.end_z === nothing
        for i in 1:nu-1
            z_ref[i+1] = z_ref[i] + slope_at(i+1) * du
        end
        return z_ref
    end
    zb = Vector{Float64}(undef, nu)
    zb[nu] = refline.end_z
    for i in nu:-1:2
        zb[i-1] = zb[i] - slope_at(i) * du
    end
    if nu >= 2
        # the true end anchor -- set directly, matching `calcRefLineZ`'s
        # `data[size-1] = info.last` at the top of the function; the blend
        # loop below deliberately never touches this index again (its bound
        # is `1:nu-2`, one short of `integrate_reference_line`'s `1:nu-1`).
        z_ref[nu] = zb[nu]
    end
    for i in 1:nu-2
        fraction = (i - 1) / (nu - 1)
        z_ref[i+1] = (1 - fraction) * (z_ref[i] + slope_at(i+1) * du) + fraction * zb[i+1]
    end
    return z_ref
end

"""
    lateral_offset_grid(x, y, v) -> (X, Y)

Offset the reference line `(x,y)` laterally by each `v` to build the full
`(nu, nv)` world-frame grid, replicating the reference implementation's
"miter normal" scheme (`crgEvaluv2xy.c`) evaluated exactly at grid nodes.

At each INTERIOR node `i`, the true C function (queried continuously between
nodes) computes an offset point by bisecting the incoming/outgoing segment
directions via the chord skipping over node `i` (from node `i-1` directly to
node `i+1`), then rescales that bisector so its projection onto the
FOLLOWING segment's own normal reproduces the true perpendicular distance
`v` (a standard "miter join" construction — this is what keeps offset
polylines from gapping/overlapping at kinks in the reference line). At
`u = u_i` exactly, only this one node's offset point is used (no
interpolation between the segment's two endpoints is needed, since
evaluating exactly at a grid node has interpolation fraction 0). The first
and last nodes have no bisector partner and fall back to their single
adjacent segment's own normal.

Perpendicular convention here: `perp(dx, dy) = (-dy, dx)` (90° CCW). If
Task 17's full-grid cross-validation shows every `v` mismatched by a sign
flip (i.e. `X`/`Y` match with `v` negated), that's this convention being
mirrored relative to the C library — fix by negating `perp`'s output here,
not by negating `v` itself, since `v`'s sign also matters for the banking
term in `assemble_z_grid` (Task 14) and must stay consistent with the
parsed `v` axis. **Confirmed during Task 17: no such flip was needed** --
cross-validating both baseline fixtures against the compiled C oracle
passed with 0 mismatches on the very first run, using this convention
exactly as written.

Epsilon-guarded degeneracy handling at an exact U-turn / coincident-endpoints
kink (found during Task 13's review; guards ported during Task 17's
cross-validation): if node `i-1` and node `i+1` coincide (or nearly so), the
"chord skipping over node `i`" has ~zero length, which would otherwise send
`normalize2` to `0/0` and the miter rescale's `denom` toward zero, producing
`NaN`/`Inf` in `offset_dir[i]` with no warning. This is NOT a purely
theoretical corner case: `integrate_reference_line`'s non-end-anchored
branch gives every segment exactly the same length `du` by construction, so
an exact-180°-hairpin reference line — plausible for a pathological/
adversarial input, if not a typical recorded track — hits this degenerate
equal-length configuration by default, not by contrived construction. The
reference implementation this function transcribes guards against exactly
this: `crgEvaluv2xy.c`'s `normalizeVector2` (~line 197) returns its input
UNCHANGED (not normalized) rather than dividing when `length < 1.0e-10`,
and its rescale step (~line 149) checks `fabs(dotProd) > 1.0e-10` before
dividing, falling back to the un-rescaled (and, per the previous guard,
possibly still un-normalized) bisector otherwise. Both thresholds are
ported verbatim below (`normalize2`/the `abs(denom) > 1.0e-10` check).
Net effect at an exact hairpin: `offset_dir[i]` ends up as a
vector with a magnitude on the order of 1e-16 (not exactly `(0.0, 0.0)`
unless the chord is bit-exactly zero-length) — so the offset point stays
essentially ON the reference line regardless of `v`, rather than blowing up.
Confirmed bit-for-bit against the compiled C oracle on a synthetic
hairpin fixture (`synthetic_hairpin.crg`, Task 17): at the hairpin node,
`crgEvaluv2xy` returns `x/y` equal to the reference line's own point there
(to ~1e-16), for every `v` queried — exactly what this guarded formula
produces. (Task 13's original pinning test asserted `isnan(X[2,1])` for
this exact scenario using EXACT `Float64` zeros for `x`/`y`, where the
chord length actually is bit-exact `0.0`, not just ~1e-16 — the guard
handles that case too, since `hypot(0.0,0.0) = 0.0 < 1.0e-10`. That test's
expectation was updated to the new, correct, finite fallback value.)
"""
function lateral_offset_grid(x::Vector{Float64}, y::Vector{Float64}, v::Vector{Float64})
    nu = length(x)
    perp(d) = (-d[2], d[1])
    # Guarded normalize: matches `normalizeVector2` (crgEvaluv2xy.c ~line 197)
    # exactly -- leave the vector UNCHANGED (not divided) when its length is
    # below 1.0e-10, instead of dividing by (near-)zero and producing NaN/Inf.
    normalize2(d) = (n = hypot(d[1], d[2]); n < 1.0e-10 ? d : (d[1]/n, d[2]/n))

    seg = [normalize2((x[i+1]-x[i], y[i+1]-y[i])) for i in 1:nu-1]
    n12 = perp.(seg)

    offset_dir = Vector{Tuple{Float64,Float64}}(undef, nu)
    if nu == 1
        offset_dir[1] = (0.0, 1.0)   # degenerate single-node "line"; arbitrary but consistent
    else
        offset_dir[1] = n12[1]
        offset_dir[nu] = n12[nu-1]
        for i in 2:nu-1
            chord = normalize2((x[i+1]-x[i-1], y[i+1]-y[i-1]))
            n1 = perp(chord)
            denom = n1[1]*n12[i][1] + n1[2]*n12[i][2]
            # Guarded rescale: matches the C reference's `fabs(dotProd) > 1.0e-10`
            # check exactly -- only rescale when the dot product is big enough;
            # otherwise fall back to the (possibly still un-normalized, per the
            # guard above) bisector `n1` as-is. This is what produces the
            # graceful "stay on the reference line" fallback at a hairpin,
            # instead of blowing up to NaN/Inf.
            offset_dir[i] = abs(denom) > 1.0e-10 ? (n1[1]/denom, n1[2]/denom) : n1
        end
    end

    nv = length(v)
    X, Y = Matrix{Float64}(undef, nu, nv), Matrix{Float64}(undef, nu, nv)
    for i in 1:nu, j in 1:nv
        X[i,j] = x[i] + v[j] * offset_dir[i][1]
        Y[i,j] = y[i] + v[j] * offset_dir[i][2]
    end
    return X, Y
end

"""
    assemble_z_grid(z_grid, z_ref, banking, refline, v) -> Z

`Z[i,j] = z_grid[i,j] + z_ref[i] + bank(i) * clamp(v[j], v[1], v[end])`.
Banking's lateral reach is clamped to the road's actual v-range regardless
of what `v` values are being queried — matching `crgEvalz.c`, which clips
`v` for the banking term specifically (but not for the grid-z lookup
itself). If there's no `banking` channel, `refline.start_banking` (default
0.0, from `REFERENCE_LINE_START_B`) is used as a constant cross-slope for
every row.

**This clamp is a provable no-op on every real call path in this package**
— UPDATED for Task 16's `apply_mods`, which now sits between
`assemble_channels` and this function in `road_surface_grid`'s actual call
chain (`d = apply_mods(data)`, then `assemble_z_grid(d.z, z_ref, d.banking,
d.refline, d.v)`). In the C reference, `crgEvaluv2z` (`crgEvalz.c` lines
~574-588) is a *continuous point-query* evaluator: its query `v` (an
arbitrary real, e.g. 5.3) is a genuinely different quantity from
`crgData->channelV.info.first`/`.last` (the channel's own declared, discrete
v-axis) — so clamping the query to the channel's declared range is real,
meaningful work. This package has no such split: it's a batched forward
transform only (see the design doc — "no per-point Newton-iteration
inversion ... forward-only"), so there is exactly one `v` in play, and
`assemble_z_grid` is always called with the same `v` its own `z_grid`
argument's columns correspond to.
That `v` traces back to `assemble_channels` (Task 9), which always returns
`v` sorted ascending (`sortperm`) with `z`'s columns permuted to match --
but by the time `road_surface_grid` calls `assemble_z_grid`, `apply_mods`
(Task 16) has had a chance to touch `v` too, via `SCALE_WIDTH`
(`v .*= mods.scale_width`). A uniform multiplicative rescale by a *positive*
factor cannot change the relative order of a sorted array, so the
ascending/matches-`z`'s-columns invariant survives (confirmed by a dedicated
regression test, `test_transform.jl`'s `"SCALE_WIDTH preserves the
v-ascending / z-column invariant..."`); this package does not guard against
a zero or negative `SCALE_WIDTH`, which would violate it — matching the C
reference, which also does not guard against that case (`crgDataScaleChannel`
in `crgMgr.c` scales `channelV`'s values unconditionally, with no sign check).
For any sorted-ascending `v`, `first(v)` and `last(v)` are exactly its min
and max, so `clamp(v[j], first(v), last(v)) == v[j]` for every `j` — always,
for any positive-`SCALE_WIDTH`-or-absent call path, not merely by
observation. The clamp is kept anyway (unconditionally, per this task's
implementation) because it's harmless, free, and mirrors the upstream
reference's defensive intent — see `test_transform.jl`'s `assemble_z_grid`
testset, which pins both (a) this no-op behavior on realistic (sorted) `v`,
and (b) the clamp's literal arithmetic via a synthetic non-monotonic `v`
that the real pipeline never produces but which still exercises the
`clamp(...)` line for regression purposes. This mirrors how Task 13's
review documented — rather than silently fixed or silently ignored — the
`lateral_offset_grid` hairpin degeneracy (since fixed in Task 17): a known,
understood limitation, pinned by a test, not swept under the rug.
"""
function assemble_z_grid(z_grid::Matrix{Float64}, z_ref::Vector{Float64}, banking::Union{Vector{Float64},Nothing}, refline::ReferenceLineParams, v::Vector{Float64})
    nu, nv = size(z_grid)
    Z = Matrix{Float64}(undef, nu, nv)
    vmin, vmax = first(v), last(v)
    for i in 1:nu
        bank_i = banking === nothing ? refline.start_banking : banking[i]
        for j in 1:nv
            vc = clamp(v[j], vmin, vmax)
            Z[i,j] = z_grid[i,j] + z_ref[i] + bank_i * vc
        end
    end
    return Z
end

"""
    apply_mods(data::CRGData) -> CRGData

Apply `\$ROAD_CRG_MODS`, returning a new `CRGData` with adjusted raw
channels. Order matches the reference implementation
(`crgDataSetModifiersApply` in `crgMgr.c`): scale channels first (z-grid,
slope, banking, length, width, curvature), then a single rotate+translate.

Key insight that keeps this simple: rotating every `phi` value by a
constant angle rotates the *shape* of the eventually-integrated reference
line by that same angle, because `(cos(a+θ), sin(a+θ)) = R(θ)·(cos a, sin a)`
— so only `refline.start_x/start_y` (rotated about the pivot, then
translated) and `phi` itself need to change here. This holds for the
`REFPOINT_*` case too: the only anchor point this function supports is the
implicit default one (`u = start_u`, `v = 0` — see the scope-boundary note
below), which is exactly `refline.start_x/start_y/start_phi` by
construction. There's no need to call `integrate_reference_line` inside
this function at all.

`REFPOINT_*` (if any of `REFPOINT_X/Y/Z/PHI` is set) takes over the
rotate+translate step entirely, ignoring `REFLINE_OFFSET_*`/
`REFLINE_ROTCENTER_*` completely — they do not compose, matching
`crgDataApplyTransformations` in `crgMgr.c`: `applyXform` is set to 1 by
any of the `dCrgModRefPoint{X,Y,Z,Phi,U,UFrac,V,VFrac}` checks, and the
whole `REFLINE_OFFSET_*`/`REFLINE_ROTCENTER_*` block is gated behind
`if (!applyXform)` — i.e. it runs at all only when NONE of those
`REFPOINT_*` fields were present. Confirmed directly against that source
(`lib/LibOpenCRG/csrc/src/crgMgr.c`, `crgDataApplyTransformations`,
~lines 804-922) rather than assumed.

**Scope boundary, narrowed during the Tasks 16+17 review (read before
assuming full `REFPOINT_*` support exists):** the `(u,v)` point a
`REFPOINT_*` modifier pins down is, in general, an arbitrary CONTINUOUS
point on the road surface — the C reference resolves it via genuine
point-interpolation (`crgDataEvaluv2xy`/`crgDataEvaluv2z` at the full query
`(uPos, vPos)`), which this package's design doc explicitly excludes from
scope (no per-point Newton-iteration inversion, no interpolation/evaluation
API, forward-only batched grid transform only). Only the ONE anchor that's
exactly computable without an interpolator is supported here: the IMPLICIT
DEFAULT point, `u = start_u` (grid node 1), `v = 0` — at which `x, y, phi`
exactly reproduce `start_x/start_y/start_phi` with no interpolation needed
(confirmed by Task 13's "v=0 reproduces the reference line exactly" test).
If `REFPOINT_U`, `REFPOINT_U_FRACTION`, `REFPOINT_V`, or
`REFPOINT_V_FRACTION` is set to ANY value at all, this function `error`s
rather than silently pinning the wrong point — matching this plan's
established "fail loudly instead of silently corrupting geometry"
philosophy (Task 9's zero-`:long_section`-channels error, Task 11's
exactly-one-of-`end_x`/`end_y` error). Relatedly, `REFPOINT_U_OFFSET`/
`REFPOINT_V_OFFSET` are only ever meaningful alongside their `_FRACTION`
sibling in the C reference (`crgMgr.c:846,865`, read only from inside the
`UFrac`/`VFrac` branches); since that whole path now always errors, they
are simply inert here, and setting either one ALONE (without its
`_FRACTION` sibling) correctly falls through to `REFLINE_OFFSET_*`/
`REFLINE_ROTCENTER_*` instead of wrongly taking the `REFPOINT_*` branch — a
real, silent, cross-validated bug in the first cut of this function (see
`test_transform.jl`'s regression test for this exact scenario).

`REFPOINT_Z`'s target elevation, like `REFPOINT_X`/`REFPOINT_Y`, is applied
even when it isn't itself set (defaulting to `0.0`, mirroring
`crgDataApplyTransformations`'s zero-initialized `toXYZ[2]`) — so setting
only `REFPOINT_X`/`REFPOINT_Y`/`REFPOINT_PHI` still shifts the reference
line's elevation profile so the anchor's elevation becomes `0.0`, unless
`REFPOINT_Z` says otherwise. The anchor's CURRENT elevation (`from_z`) is
NOT just `z_ref[1]` (the reference-line elevation channel alone) —
cross-validating a `REFPOINT_Z`-only scenario against the real C oracle
(`test_transform.jl`'s `"REFPOINT_Z"` testset, using `belgian_block.crg`,
whose z-grid has a real ~2.13m nonzero value at `(start_u, v=0)`) showed
`from_z` must fold in the z-grid's own value at `v = 0` too, i.e.
`from_z = z_ref[1] + z_grid[1, v0_idx]`, matching `crgDataEvaluv2z`'s full
surface evaluation (grid + reference elevation + banking, the last of which
vanishes at `v=0`) rather than the reference-elevation-only value. This
requires an EXACT `v = 0` entry in the (possibly `SCALE_WIDTH`-adjusted)
`v` axis to look up without interpolating; if none exists, this function
`error`s rather than approximating (same philosophy as the
`REFPOINT_U*`/`REFPOINT_V*` error above).

**`SCALE_LENGTH` also rescales `end_u` proportionally about `start_u`, not just
`increment`** — matching `crgDataSetModifiersApply`'s explicit
`channelU.info.last = channelU.info.first + dValue * uRange` (`crgMgr.c`
~lines 548-558), where `uRange` is the ORIGINAL (pre-scale) `end_u - start_u`,
and `channelU.info.first` (`start_u`) is left untouched. Found during the
whole-plan final review: an earlier version of this function scaled
`u_increment` correctly but threaded `r.end_u` straight through to `r2`/`r3`
completely unchanged. `road_surface_grid`'s own u-axis (`u = [start_u +
(i-1)*increment for i in 1:nu]`) only ever reads `start_u`/`increment`, never
`end_u` — so this had (and, now fixed, still has) zero effect on any grid
this package actually produces. But `end_u` IS part of the public
`ReferenceLineParams`/`CRGData` struct a caller can read directly, so leaving
it stale after a length rescale was a real correctness bug regardless of
downstream reachability — cross-validated directly against the compiled C
oracle's own post-scale `crgDataSetGetURange` output (not just hand
arithmetic) in `test_transform.jl`'s `"SCALE_LENGTH"` regression test.

Deliberately NOT touched by any other `SCALE_*` here: `v` only moves under
`SCALE_WIDTH`, and only by a uniform multiplicative factor — this preserves
(assuming a positive scale factor, the only physically sensible one for a
road width) the invariant `assemble_z_grid` (Task 14) depends on, that `v`
stays sorted ascending and still lines up 1:1 with `z`'s columns (Task 9's
`sortperm` guarantee is about relative order, which scaling by a positive
constant can't disturb). See `test_transform.jl`'s
`"SCALE_WIDTH preserves the v-ascending / z-column invariant..."` regression
test.

Deliberately NOT applied: `mods.grid_nan_mode`/`mods.grid_nan_offset` (parsed
by Task 6, real behavior in `crgMgr.c` around line 601) control how NaN gaps
in the z-grid get filled/replaced — a data-cleaning concern, not the
scale/offset/rotate geometry transform this task's scope was explicitly
limited to. Same deferral treatment as `\$ROAD_CRG_MPRO` (Task 5's design
doc note) — parsed and available on `RoadCrgMods`, never auto-applied.
"""
function apply_mods(data::CRGData)
    mods = data.mods
    r = data.refline
    phi = copy(data.phi)
    z = copy(data.z)
    slope = data.slope === nothing ? nothing : copy(data.slope)
    banking = data.banking === nothing ? nothing : copy(data.banking)
    v = copy(data.v)
    u_increment, start_slope, end_slope, start_banking, end_banking =
        r.increment, r.start_slope, r.end_slope, r.start_banking, r.end_banking

    mods.scale_z_grid === nothing || (z .*= mods.scale_z_grid)
    if mods.scale_slope !== nothing
        slope === nothing || (slope .*= mods.scale_slope)
        start_slope *= mods.scale_slope
        end_slope = end_slope === nothing ? nothing : end_slope * mods.scale_slope
    end
    if mods.scale_banking !== nothing
        banking === nothing || (banking .*= mods.scale_banking)
        start_banking *= mods.scale_banking
        end_banking = end_banking === nothing ? nothing : end_banking * mods.scale_banking
    end
    mods.scale_length === nothing || (u_increment *= mods.scale_length)
    mods.scale_width === nothing || (v .*= mods.scale_width)
    if mods.scale_curvature !== nothing
        base = phi[1]
        for i in 2:length(phi)
            phi[i] = base + mods.scale_curvature * (phi[i] - base)
        end
    end

    r2 = ReferenceLineParams(r.start_u, r.end_u, u_increment, r.start_x, r.start_y, r.start_phi,
        r.end_x, r.end_y, r.end_phi, r.start_z, r.end_z, r.v_right, r.v_left, r.v_increment,
        start_slope, end_slope, start_banking, end_banking)

    # SCALE_LENGTH: end_u must rescale proportionally about start_u (matching
    # crgDataSetModifiersApply's `channelU.info.last = channelU.info.first +
    # dValue * uRange`, crgMgr.c ~548-558) -- NOT stay fixed at r2.end_u, which
    # is otherwise never read again between here and r3's construction below
    # (only threaded straight through), making this the one place the
    # corrected value needs to land. See this function's docstring.
    new_end_u = r2.start_u + something(mods.scale_length, 1.0) * (r2.end_u - r2.start_u)

    # REFPOINT_U/REFPOINT_U_FRACTION/REFPOINT_V/REFPOINT_V_FRACTION pin an
    # arbitrary CONTINUOUS (u,v) point, which requires genuine point
    # interpolation (crgDataEvaluv2xy/crgDataEvaluv2z in the C reference) --
    # out of scope for this package's grid-only, forward-only design (see
    # this function's docstring). Error unconditionally on ANY value, rather
    # than special-casing "happens to equal the default", which would be
    # fragile, floating-point-exactness-dependent behavior.
    any(f -> getfield(mods, f) !== nothing,
        (:refpoint_u, :refpoint_u_fraction, :refpoint_v, :refpoint_v_fraction)) &&
        error("apply_mods: REFPOINT_U/REFPOINT_U_FRACTION/REFPOINT_V/REFPOINT_V_FRACTION are not " *
              "supported -- pinning an arbitrary continuous (u,v) point requires per-point " *
              "interpolation, which this package's grid-only, forward-only design does not " *
              "implement. Only REFPOINT_X/REFPOINT_Y/REFPOINT_Z/REFPOINT_PHI, anchored at the " *
              "implicit default point (u = start_u, v = 0), are supported.")

    # NOTE: REFPOINT_U_OFFSET/REFPOINT_V_OFFSET are deliberately absent from
    # this trigger tuple -- in the C reference they are only ever read from
    # INSIDE the UFrac/VFrac branches (crgMgr.c:846,865), never independent
    # triggers, so (matching that exactly) they must be inert here too, now
    # that the UFrac/VFrac path they'd otherwise feed always errors above.
    has_refpoint = any(f -> getfield(mods, f) !== nothing,
        (:refpoint_x, :refpoint_y, :refpoint_z, :refpoint_phi))

    if has_refpoint
        # The only supported anchor is the implicit default (u = start_u,
        # i.e. grid node 1; v = 0), which is exactly the reference line's own
        # start point -- no integration/interpolation needed for x/y/phi.
        from_x, from_y, from_phi = r2.start_x, r2.start_y, r2.start_phi
        rot_center = (from_x, from_y)
        rot_angle = something(mods.refpoint_phi, 0.0) - from_phi
        translation = (something(mods.refpoint_x, 0.0) - from_x, something(mods.refpoint_y, 0.0) - from_y)

        # Matching crgDataApplyTransformations exactly: the anchor's current
        # elevation (fromXYZ[2] there) is evaluated via crgDataEvaluv2z's
        # FULL z surface -- z-grid value at v=0 PLUS the reference-line
        # elevation, banking excluded since it's multiplied by v=0 -- not the
        # reference-line elevation alone, and this happens whenever ANY of
        # REFPOINT_X/Y/Z/PHI triggers this branch, not only when REFPOINT_Z
        # itself is set (cross-validated; see the docstring and
        # test_transform.jl's "REFPOINT_Z" testset).
        v0_idx = findfirst(==(0.0), v)
        v0_idx === nothing && error("apply_mods: REFPOINT_X/Y/Z/PHI need the elevation at the default " *
              "anchor's v=0, which requires an exact v=0 declared long-section column -- this package " *
              "does not implement v-interpolation, and this file has no exact v=0 column.")
        # The reference-line component of the anchor's elevation is just
        # `r2.start_z`, not `integrate_reference_z(r2, slope, length(phi))[1]`
        # -- BOTH of integrate_reference_z's paths (the early-return AND the
        # general path) set element 1 to `refline.start_z` unconditionally
        # before any loop/blend logic runs (see its implementation), so that
        # element is always exactly `r2.start_z`, regardless of slope/end_z.
        # Calling the full O(nu) integration (including its forward/backward
        # + blend double loop when end_z is set) just to read element 1 is
        # wasted work; reading `r2.start_z` directly is identical.
        from_z = r2.start_z + z[1, v0_idx]
        z_shift = something(mods.refpoint_z, 0.0) - from_z
    else
        rot_center = (something(mods.refline_rotcenter_x, r2.start_x), something(mods.refline_rotcenter_y, r2.start_y))
        rot_angle = something(mods.refline_offset_phi, 0.0)
        translation = (something(mods.refline_offset_x, 0.0), something(mods.refline_offset_y, 0.0))
        z_shift = something(mods.refline_offset_z, 0.0)
    end

    c, s = cos(rot_angle), sin(rot_angle)
    dx, dy = r2.start_x - rot_center[1], r2.start_y - rot_center[2]
    new_start_x = rot_center[1] + dx*c - dy*s + translation[1]
    new_start_y = rot_center[2] + dx*s + dy*c + translation[2]
    new_start_phi = r2.start_phi + rot_angle
    phi .+= rot_angle
    new_start_z = r2.start_z + z_shift
    new_end_z = r2.end_z === nothing ? nothing : r2.end_z + z_shift
    new_end_x = r2.end_x === nothing ? nothing : rot_center[1] + (r2.end_x-rot_center[1])*c - (r2.end_y-rot_center[2])*s + translation[1]
    new_end_y = r2.end_y === nothing ? nothing : rot_center[2] + (r2.end_x-rot_center[1])*s + (r2.end_y-rot_center[2])*c + translation[2]
    new_end_phi = r2.end_phi === nothing ? nothing : r2.end_phi + rot_angle

    r3 = ReferenceLineParams(r2.start_u, new_end_u, r2.increment, new_start_x, new_start_y, new_start_phi,
        new_end_x, new_end_y, new_end_phi, new_start_z, new_end_z,
        r2.v_right, r2.v_left, r2.v_increment, r2.start_slope, r2.end_slope, r2.start_banking, r2.end_banking)

    return CRGData(data.comment, r3, data.format_code, data.opts, mods, data.mpro, phi, banking, slope, v, z)
end

"""
    road_surface_grid(data::CRGData) -> (u, v, X, Y, Z)

The batched forward transform: applies `\$ROAD_CRG_MODS` (a no-op if none
are set), integrates the reference line and its elevation profile, and
offsets the whole grid laterally — producing world-frame `(u, v, X, Y, Z)`,
all `(nu, nv)` except the `u`/`v` axis vectors themselves.
"""
function road_surface_grid(data::CRGData)
    d = apply_mods(data)
    nu = length(d.phi)
    u = [d.refline.start_u + (i-1)*d.refline.increment for i in 1:nu]
    x, y = integrate_reference_line(d.refline, d.phi)
    z_ref = integrate_reference_z(d.refline, d.slope, nu)
    X, Y = lateral_offset_grid(x, y, d.v)
    Z = assemble_z_grid(d.z, z_ref, d.banking, d.refline, d.v)
    return u, d.v, X, Y, Z
end
