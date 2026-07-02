using OpenCRG

const _BELGIAN_BLOCK_GRID_CACHE = Ref{Any}(nothing)

function belgian_block_crg_path()
    return joinpath(dirname(@__DIR__), "lib", "OpenCRG", "test", "data", "belgian_block.crg")
end

function belgian_block_road_surface_grid()
    cached = _BELGIAN_BLOCK_GRID_CACHE[]
    if cached !== nothing
        return cached
    end
    data = OpenCRG.read_crg(belgian_block_crg_path())
    grid = OpenCRG.road_surface_grid(data)
    _BELGIAN_BLOCK_GRID_CACHE[] = grid
    return grid
end

function _padded_axis_and_values(axis::AbstractVector{<:Real}, values::AbstractVector{<:Real}; padding::Real = 5.0)
    @assert length(axis) == length(values) "Axis and value vectors must have the same length."
    @assert length(axis) >= 2 "Need at least two points to infer road-profile spacing."
    @assert all(>(0), diff(axis)) "Road-profile axis must be strictly increasing."
    @assert all(isfinite, values) "Road-profile values must be finite."

    axis_f = collect(Float64, axis)
    values_f = collect(Float64, values)
    dx = minimum(diff(axis_f))
    n_pad = max(1, ceil(Int, padding / dx))
    left_axis = collect((first(axis_f) - n_pad * dx):dx:(first(axis_f) - dx))
    right_axis = collect((last(axis_f) + dx):dx:(last(axis_f) + n_pad * dx))
    return vcat(left_axis, axis_f, right_axis), vcat(fill(first(values_f), length(left_axis)), values_f, fill(last(values_f), length(right_axis)))
end

"""
    belgian_block_centerline_profile(; degree = 3, padding = 5.0)

Return `(x, z)` for the Belgian-block centerline road profile, shifted so the
longitudinal coordinate starts at zero and the first centerline height is zero.
The profile is padded with constant end values so a vehicle can start slightly
before the measured CRG segment.
"""
function belgian_block_centerline_profile(; padding::Real = 5.0)
    u, v, _X, _Y, Z = belgian_block_road_surface_grid()
    center_j = argmin(abs.(v .- 0.0))
    z = collect(Float64, Z[:, center_j])
    finite_mask = isfinite.(z)
    x = collect(Float64, u .- first(u))
    x = x[finite_mask]
    z = z[finite_mask]
    z .-= first(z)
    return _padded_axis_and_values(x, z; padding)
end

"""
    belgian_block_centerline_profile_interpolator(; degree = 3, padding = 5.0, max_derivative_order_eval = 1)

Build a smooth one-dimensional B-spline callable for the Belgian-block centerline
height profile. The returned function maps local longitudinal position in meters
to relative road height in meters.
"""
function belgian_block_centerline_profile_interpolator(; degree::Integer = 3, padding::Real = 5.0, max_derivative_order_eval::Integer = 1)
    x, z = belgian_block_centerline_profile(; padding)
    return nd_bspline_interpolation((x,), z; degree = degree, max_derivative_order_eval = max_derivative_order_eval)
end

belgian_block_centerline_profile_value(x) = belgian_block_centerline_profile_interpolator()(x)
