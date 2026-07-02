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
