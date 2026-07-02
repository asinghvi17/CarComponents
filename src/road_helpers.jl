using OpenCRG

function belgian_block_crg_path()
    return joinpath(dirname(@__DIR__), "lib", "OpenCRG", "test", "data", "belgian_block.crg")
end

function opencrg_road_surface_grid(source)
    if source isa AbstractString
        return OpenCRG.road_surface_grid(OpenCRG.read_crg(abspath(source)))
    elseif source isa OpenCRG.CRGData
        return OpenCRG.road_surface_grid(source)
    elseif source isa Tuple && length(source) == 5
        return source
    else
        error("Expected an OpenCRG file path, OpenCRG.CRGData, or a (u, v, X, Y, Z) grid tuple.")
    end
end

belgian_block_road_surface_grid() = opencrg_road_surface_grid(belgian_block_crg_path())

function _finite_lateral_columns(Z)
    finite_cols = [all(isfinite, Z[:, j]) for j in axes(Z, 2)]
    first_col = findfirst(finite_cols)
    last_col = findlast(finite_cols)
    @assert first_col !== nothing && last_col !== nothing "No finite lateral road-surface band found."
    return first_col:last_col
end

function _select_lateral_columns(v, finite_cols; lateral_min = nothing, lateral_max = nothing)
    cols = collect(finite_cols)
    if lateral_min !== nothing
        cols = filter(j -> v[j] >= lateral_min, cols)
    end
    if lateral_max !== nothing
        cols = filter(j -> v[j] <= lateral_max, cols)
    end
    @assert !isempty(cols) "No finite lateral columns remain after applying lateral limits."
    return first(cols):last(cols)
end

"""
    opencrg_centerline_profile(source; lateral = 0.0)

Return `(x, z)` for a one-dimensional road-height profile from any OpenCRG
source. `source` can be a CRG file path, `OpenCRG.CRGData`, or a precomputed
`(u, v, X, Y, Z)` grid tuple. The longitudinal coordinate is shifted so the CRG
start is `x = 0`, and height is shifted so the first selected point is zero.
"""
function opencrg_centerline_profile(source; lateral::Real = 0.0)
    u, v, _X, _Y, Z = opencrg_road_surface_grid(source)
    j = argmin(abs.(v .- lateral))
    heights = collect(Float64, Z[:, j])
    finite_mask = isfinite.(heights)
    x = collect(Float64, u .- first(u))
    x = x[finite_mask]
    heights = heights[finite_mask]
    heights .-= first(heights)
    return x, heights
end

"""
    opencrg_centerline_profile_interpolator(source; lateral = 0.0, degree = 3, max_derivative_order_eval = 1)

Build a one-dimensional B-spline road-height profile from any OpenCRG source.
Out-of-range lookup behavior follows DataInterpolationsND's B-spline boundary
behavior.
"""
function opencrg_centerline_profile_interpolator(source; lateral::Real = 0.0, degree::Integer = 3, max_derivative_order_eval::Integer = 1)
    x, heights = opencrg_centerline_profile(source; lateral)
    return nd_bspline_interpolation((x,), heights; degree, max_derivative_order_eval)
end

belgian_block_centerline_profile(; kwargs...) = opencrg_centerline_profile(belgian_block_crg_path(); kwargs...)
belgian_block_centerline_profile_interpolator(; kwargs...) = opencrg_centerline_profile_interpolator(belgian_block_crg_path(); kwargs...)
belgian_block_centerline_profile_value(x) = belgian_block_centerline_profile_interpolator()(x)

"""
    opencrg_surface_data(source; lateral_min = nothing, lateral_max = nothing, reference_lateral = 0.0)

Return `(x_axis, z_axis, height)` for a finite OpenCRG road surface in the car's
world X-Z driving plane. X is local longitudinal distance from the CRG start, Z
is the road lateral coordinate, and height is relative to the road height at the
first longitudinal station and nearest `reference_lateral`.

Rows/columns are not padded. Out-of-range lookup behavior is left to
DataInterpolationsND.
"""
function opencrg_surface_data(source; lateral_min = nothing, lateral_max = nothing, reference_lateral::Real = 0.0)
    u, v, _X, _Y, Z = opencrg_road_surface_grid(source)
    finite_cols = _finite_lateral_columns(Z)
    cols = _select_lateral_columns(v, finite_cols; lateral_min, lateral_max)

    x_axis = collect(Float64, u .- first(u))
    z_axis = collect(Float64, v[cols])
    heights = collect(Float64, Z[:, cols])
    @assert all(isfinite, heights) "Selected OpenCRG surface contains non-finite height values."

    reference_col = argmin(abs.(z_axis .- reference_lateral))
    heights .-= heights[1, reference_col]

    return x_axis, z_axis, heights
end

"""
    opencrg_surface_interpolator(source; degree_x = 1, degree_z = 1, max_derivative_order_eval = 1, ...)

Build a 2D B-spline road-height surface from any OpenCRG source. The returned
callable maps `(world_x, world_z)` to road height in world Y, using the car's X-Z
driving plane. Extrapolation follows DataInterpolationsND's B-spline boundary
behavior.
"""
function opencrg_surface_interpolator(
        source;
        degree_x::Integer = 1,
        degree_z::Integer = 1,
        max_derivative_order_eval::Integer = 1,
        lateral_min = nothing,
        lateral_max = nothing,
        reference_lateral::Real = 0.0,
    )
    x_axis, z_axis, heights = opencrg_surface_data(source; lateral_min, lateral_max, reference_lateral)
    return nd_bspline_interpolation(
        (x_axis, z_axis),
        heights;
        degree = (degree_x, degree_z),
        max_derivative_order_eval = max_derivative_order_eval,
    )
end

belgian_block_surface_data(; kwargs...) = opencrg_surface_data(belgian_block_crg_path(); kwargs...)
belgian_block_surface_interpolator(; kwargs...) = opencrg_surface_interpolator(belgian_block_crg_path(); kwargs...)
belgian_block_surface_value(x, z) = belgian_block_surface_interpolator()(x, z)
