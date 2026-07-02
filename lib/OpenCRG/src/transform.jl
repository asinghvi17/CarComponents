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
    for i in 1:nu-1
        fraction = i / (nu - 1)
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
parsed `v` axis.

Known unguarded degeneracy (found during Task 13's review, NOT fixed here —
deferred to Task 17): an exact U-turn / coincident-endpoints kink (node
`i-1` and node `i+1` at the same point, so the "chord skipping over node i"
has zero length) sends `normalize2` to `0/0` and the miter rescale's
`denom` toward zero, producing `NaN`/`Inf` in `offset_dir[i]` with no
warning. This is NOT a purely theoretical corner case: `integrate_reference_line`'s
non-end-anchored branch gives every segment exactly the same length `du` by
construction, so an exact-180°-hairpin reference line — plausible for a
pathological/adversarial input, if not a typical recorded track — hits this
degenerate equal-length configuration by default, not by contrived
construction. The reference implementation this task transcribes already
guards against exactly this: `crgEvaluv2xy.c`'s `normalizeVector2` (~line
197) returns its input unchanged rather than dividing when
`length < 1.0e-10`, and its rescale step (~line 149) checks
`fabs(dotProd) > 1.0e-10` before dividing, falling back to the un-rescaled
normal otherwise. Task 17 should port those exact thresholds/fallbacks
(cross-validatable directly against the compiled C oracle via FFI) rather
than deriving new ones from scratch.
"""
function lateral_offset_grid(x::Vector{Float64}, y::Vector{Float64}, v::Vector{Float64})
    nu = length(x)
    perp(d) = (-d[2], d[1])
    normalize2(d) = (n = hypot(d[1], d[2]); (d[1]/n, d[2]/n))

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
            offset_dir[i] = (n1[1]/denom, n1[2]/denom)
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
translated) and `phi` itself need to change here. There's no need to
integrate the reference line inside this function at all, EXCEPT for the
`REFPOINT_*` case, which needs to evaluate the CURRENT (pre-transform)
position at a specific `(u,v)` to know what point is being pinned down.

`REFPOINT_*` (if ANY such field is set) takes over the rotate+translate
step entirely, ignoring `REFLINE_OFFSET_*`/`REFLINE_ROTCENTER_*` completely
— they do not compose, matching `crgDataApplyTransformations` in `crgMgr.c`:
`applyXform` is set to 1 by any of the `dCrgModRefPoint{X,Y,Z,Phi,U,UFrac,V,VFrac}`
checks, and the whole `REFLINE_OFFSET_*`/`REFLINE_ROTCENTER_*` block is
gated behind `if (!applyXform)` — i.e. it runs at all only when NONE of
those `REFPOINT_*` fields were present. Confirmed directly against that
source (`lib/LibOpenCRG/csrc/src/crgMgr.c`, `crgDataApplyTransformations`,
~lines 804-922) rather than assumed.

Deliberately NOT touched by any `SCALE_*` here: `v` only moves under
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

    has_refpoint = any(f -> getfield(mods, f) !== nothing, (:refpoint_u, :refpoint_u_fraction,
        :refpoint_u_offset, :refpoint_v, :refpoint_v_fraction, :refpoint_v_offset,
        :refpoint_x, :refpoint_y, :refpoint_z, :refpoint_phi))

    if has_refpoint
        x0, y0 = integrate_reference_line(r2, phi)
        u_frac, u_off = something(mods.refpoint_u_fraction, 0.0), something(mods.refpoint_u_offset, 0.0)
        u_pos = something(mods.refpoint_u, r2.start_u + u_frac * (r2.end_u - r2.start_u) + u_off)
        idx = clamp(round(Int, (u_pos - r2.start_u) / u_increment) + 1, 1, length(phi))
        from_x, from_y, from_phi = x0[idx], y0[idx], phi[idx]
        rot_center = (from_x, from_y)
        rot_angle = something(mods.refpoint_phi, 0.0) - from_phi
        translation = (something(mods.refpoint_x, 0.0) - from_x, something(mods.refpoint_y, 0.0) - from_y)
    else
        rot_center = (something(mods.refline_rotcenter_x, r2.start_x), something(mods.refline_rotcenter_y, r2.start_y))
        rot_angle = something(mods.refline_offset_phi, 0.0)
        translation = (something(mods.refline_offset_x, 0.0), something(mods.refline_offset_y, 0.0))
    end

    c, s = cos(rot_angle), sin(rot_angle)
    dx, dy = r2.start_x - rot_center[1], r2.start_y - rot_center[2]
    new_start_x = rot_center[1] + dx*c - dy*s + translation[1]
    new_start_y = rot_center[2] + dx*s + dy*c + translation[2]
    new_start_phi = r2.start_phi + rot_angle
    phi .+= rot_angle
    new_start_z = r2.start_z + something(mods.refline_offset_z, 0.0)
    new_end_z = r2.end_z === nothing ? nothing : r2.end_z + something(mods.refline_offset_z, 0.0)
    new_end_x = r2.end_x === nothing ? nothing : rot_center[1] + (r2.end_x-rot_center[1])*c - (r2.end_y-rot_center[2])*s + translation[1]
    new_end_y = r2.end_y === nothing ? nothing : rot_center[2] + (r2.end_x-rot_center[1])*s + (r2.end_y-rot_center[2])*c + translation[2]
    new_end_phi = r2.end_phi === nothing ? nothing : r2.end_phi + rot_angle

    r3 = ReferenceLineParams(r2.start_u, r2.end_u, r2.increment, new_start_x, new_start_y, new_start_phi,
        new_end_x, new_end_y, new_end_phi, new_start_z, new_end_z,
        r2.v_right, r2.v_left, r2.v_increment, r2.start_slope, r2.end_slope, r2.start_banking, r2.end_banking)

    return CRGData(data.comment, r3, data.format_code, data.opts, mods, data.mpro, phi, banking, slope, v, z)
end
