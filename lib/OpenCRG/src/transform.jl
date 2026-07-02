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
