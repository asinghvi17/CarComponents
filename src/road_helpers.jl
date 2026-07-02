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

function _pad_surface_axes_and_values(
        x_axis::AbstractVector{<:Real},
        z_axis::AbstractVector{<:Real},
        values::AbstractMatrix{<:Real};
        x_padding::Real = 5.0,
        z_padding::Real = 0.5,
    )
    @assert size(values) == (length(x_axis), length(z_axis)) "Surface matrix size must match axis lengths."
    @assert all(>(0), diff(x_axis)) "Surface x-axis must be strictly increasing."
    @assert all(>(0), diff(z_axis)) "Surface z-axis must be strictly increasing."
    @assert all(isfinite, values) "Surface values must be finite."

    x = collect(Float64, x_axis)
    z = collect(Float64, z_axis)
    Z = collect(Float64, values)

    dx = minimum(diff(x))
    nx_left = max(1, ceil(Int, x_padding / dx))
    nx_right = nx_left
    x_left = collect((first(x) - nx_left * dx):dx:(first(x) - dx))
    x_right = collect((last(x) + dx):dx:(last(x) + nx_right * dx))
    x_pad = vcat(x_left, x, x_right)
    Zx = vcat(repeat(Z[1:1, :], nx_left, 1), Z, repeat(Z[end:end, :], nx_right, 1))

    dz = minimum(diff(z))
    nz_left = max(1, ceil(Int, z_padding / dz))
    nz_right = nz_left
    z_left = collect((first(z) - nz_left * dz):dz:(first(z) - dz))
    z_right = collect((last(z) + dz):dz:(last(z) + nz_right * dz))
    z_pad = vcat(z_left, z, z_right)
    Zpad = hcat(repeat(Zx[:, 1:1], 1, nz_left), Zx, repeat(Zx[:, end:end], 1, nz_right))

    return x_pad, z_pad, Zpad
end

"""
    belgian_block_centerline_profile(; padding = 5.0)

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

"""
    belgian_block_surface_data(; x_padding = 5.0, z_padding = 0.5)

Return `(x_axis, z_axis, height)` for a finite Belgian-block road surface in the
car's world X-Z driving plane. X is the local longitudinal distance from the CRG
start, Z is the road lateral coordinate, and height is relative to the centerline
height at the CRG start.
"""
function belgian_block_surface_data(; x_padding::Real = 5.0, z_padding::Real = 0.5)
    u, v, _X, _Y, Z = belgian_block_road_surface_grid()
    finite_cols = [all(isfinite, Z[:, j]) for j in axes(Z, 2)]
    first_col = findfirst(finite_cols)
    last_col = findlast(finite_cols)
    @assert first_col !== nothing && last_col !== nothing "No finite lateral road-surface band found."

    x_axis = collect(Float64, u .- first(u))
    z_axis = collect(Float64, v[first_col:last_col])
    heights = collect(Float64, Z[:, first_col:last_col])

    center_j = argmin(abs.(z_axis .- 0.0))
    heights .-= heights[1, center_j]

    return _pad_surface_axes_and_values(x_axis, z_axis, heights; x_padding, z_padding)
end

"""
    belgian_block_surface_interpolator(; degree_x = 3, degree_z = 3, x_padding = 5.0, z_padding = 0.5)

Build a smooth 2D B-spline road-height surface. The returned callable maps
`(world_x, world_z)` to road height in world Y, using the car's X-Z driving plane.
"""
function belgian_block_surface_interpolator(; degree_x::Integer = 3, degree_z::Integer = 3, x_padding::Real = 5.0, z_padding::Real = 0.5, max_derivative_order_eval::Integer = 1)
    x_axis, z_axis, heights = belgian_block_surface_data(; x_padding, z_padding)
    return nd_bspline_interpolation(
        (x_axis, z_axis),
        heights;
        degree = (degree_x, degree_z),
        max_derivative_order_eval = max_derivative_order_eval,
    )
end

function _finite_difference_along_axis(values::AbstractMatrix{<:Real}, axis::AbstractVector{<:Real}, dim::Integer)
    slopes = similar(collect(Float64, values))
    if dim == 1
        for j in axes(values, 2), i in axes(values, 1)
            if i == firstindex(axis)
                slopes[i, j] = (values[i + 1, j] - values[i, j]) / (axis[i + 1] - axis[i])
            elseif i == lastindex(axis)
                slopes[i, j] = (values[i, j] - values[i - 1, j]) / (axis[i] - axis[i - 1])
            else
                slopes[i, j] = (values[i + 1, j] - values[i - 1, j]) / (axis[i + 1] - axis[i - 1])
            end
        end
    elseif dim == 2
        for i in axes(values, 1), j in axes(values, 2)
            if j == firstindex(axis)
                slopes[i, j] = (values[i, j + 1] - values[i, j]) / (axis[j + 1] - axis[j])
            elseif j == lastindex(axis)
                slopes[i, j] = (values[i, j] - values[i, j - 1]) / (axis[j] - axis[j - 1])
            else
                slopes[i, j] = (values[i, j + 1] - values[i, j - 1]) / (axis[j + 1] - axis[j - 1])
            end
        end
    else
        error("dim must be 1 or 2")
    end
    return slopes
end

function belgian_block_surface_slope_x_interpolator(; degree_x::Integer = 3, degree_z::Integer = 3, x_padding::Real = 5.0, z_padding::Real = 0.5, max_derivative_order_eval::Integer = 1)
    x_axis, z_axis, heights = belgian_block_surface_data(; x_padding, z_padding)
    slope_x = _finite_difference_along_axis(heights, x_axis, 1)
    return nd_bspline_interpolation(
        (x_axis, z_axis),
        slope_x;
        degree = (degree_x, degree_z),
        max_derivative_order_eval = max_derivative_order_eval,
    )
end

function belgian_block_surface_slope_z_interpolator(; degree_x::Integer = 3, degree_z::Integer = 3, x_padding::Real = 5.0, z_padding::Real = 0.5, max_derivative_order_eval::Integer = 1)
    x_axis, z_axis, heights = belgian_block_surface_data(; x_padding, z_padding)
    slope_z = _finite_difference_along_axis(heights, z_axis, 2)
    return nd_bspline_interpolation(
        (x_axis, z_axis),
        slope_z;
        degree = (degree_x, degree_z),
        max_derivative_order_eval = max_derivative_order_eval,
    )
end

belgian_block_surface_value(x, z) = belgian_block_surface_interpolator()(x, z)

"""
    surface_orientation_from_slopes(slope_x, slope_z)

Rotation matrix for a local road tangent plane `y = f(x, z)` using precomputed
slopes. It maps `surface_frame` axes so +Y is the local road normal and +X is the
road longitudinal tangent. Flat road gives the identity matrix.
"""
function surface_orientation_from_slopes(slope_x, slope_z)
    a = sqrt(1 + slope_x * slope_x)
    n = sqrt(1 + slope_x * slope_x + slope_z * slope_z)
    return [
        1 / a                         slope_x / a   0
        -slope_x / n                  1 / n        -slope_z / n
        -slope_x * slope_z / (a * n)  slope_z / (a * n)  a / n
    ]
end
