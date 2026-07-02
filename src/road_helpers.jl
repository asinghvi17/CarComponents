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


function country_road_crg_path()
    return joinpath(dirname(@__DIR__), "scripts", "data", "real_roads", "country_road.crg")
end

country_road_surface_grid() = opencrg_road_surface_grid(country_road_crg_path())

function _opencrg_row_stride(u, longitudinal_spacing)
    longitudinal_spacing === nothing && return 1
    du = (last(u) - first(u)) / max(length(u) - 1, 1)
    return max(1, round(Int, longitudinal_spacing / du))
end

function _opencrg_lateral_stride(lateral_stride)
    lateral_stride === nothing && return 1
    return max(1, Int(lateral_stride))
end

function _include_last_index(idx, last_index::Integer)
    values = collect(idx)
    if values[end] != last_index
        push!(values, last_index)
    end
    return values
end

function _transformed_opencrg_xy(X, Y, reference_col)
    x0 = X[1, reference_col]
    y0 = Y[1, reference_col]
    heading0 = atan(Y[2, reference_col] - Y[1, reference_col], X[2, reference_col] - X[1, reference_col])
    c = cos(heading0)
    s = sin(heading0)

    Xlocal = similar(X, Float64)
    Zlocal = similar(Y, Float64)
    for I in eachindex(X)
        dx = X[I] - x0
        dy = Y[I] - y0
        Xlocal[I] = c * dx + s * dy
        Zlocal[I] = -s * dx + c * dy
    end

    return Xlocal, Zlocal, heading0
end

function _centerline_heading_from_xz(x, z)
    heading = similar(x, Float64)
    for i in eachindex(x)
        im = max(firstindex(x), i - 1)
        ip = min(lastindex(x), i + 1)
        heading[i] = atan(z[ip] - z[im], x[ip] - x[im])
    end
    return heading
end

"""
    opencrg_curved_road_data(source; longitudinal_spacing = 0.25, lateral_stride = 4, reference_lateral = 0.0)

Return a downsampled, transformed OpenCRG road data set for simulation and
animation. The returned X/Z coordinates are world path-plane coordinates with the
reference line start translated to `(0, 0)` and rotated so the initial heading is
along +X. Height is shifted relative to the first reference-line sample.
"""
function opencrg_curved_road_data(
        source;
        longitudinal_spacing = 0.25,
        lateral_stride = 4,
        lateral_min = nothing,
        lateral_max = nothing,
        reference_lateral::Real = 0.0,
    )
    u, v, X, Y, Z = opencrg_road_surface_grid(source)
    finite_cols = _finite_lateral_columns(Z)
    selected_cols = _select_lateral_columns(v, finite_cols; lateral_min, lateral_max)
    reference_col = argmin(abs.(v .- reference_lateral))
    @assert reference_col in selected_cols "The reference lateral column must be finite and inside the selected lateral band."

    row_stride = _opencrg_row_stride(u, longitudinal_spacing)
    col_stride = _opencrg_lateral_stride(lateral_stride)
    rows = _include_last_index(1:row_stride:length(u), length(u))
    cols = _include_last_index(first(selected_cols):col_stride:last(selected_cols), last(selected_cols))
    if !(reference_col in cols)
        push!(cols, reference_col)
        sort!(cols)
    end

    Xlocal_full, Zlocal_full, initial_heading = _transformed_opencrg_xy(X, Y, reference_col)
    x_axis = collect(Float64, Xlocal_full[rows, reference_col])
    center_z = collect(Float64, Zlocal_full[rows, reference_col])
    center_heading = _centerline_heading_from_xz(x_axis, center_z)
    v_axis = collect(Float64, v[cols])
    reference_height = Z[1, reference_col]
    heights = collect(Float64, Z[rows, cols] .- reference_height)
    @assert all(isfinite, heights) "Selected OpenCRG road data contains non-finite height values."
    @assert all(>(0), diff(x_axis)) "Transformed reference-line X coordinate must be strictly increasing."

    return (
        x_axis = x_axis,
        v_axis = v_axis,
        heights = heights,
        center_z = center_z,
        center_heading = center_heading,
        X = collect(Float64, Xlocal_full[rows, cols]),
        Z = collect(Float64, Zlocal_full[rows, cols]),
        Y = heights,
        u = collect(Float64, u[rows] .- first(u)),
        initial_heading = initial_heading,
        reference_height = reference_height,
    )
end

country_road_curved_data(; kwargs...) = opencrg_curved_road_data(country_road_crg_path(); kwargs...)

function opencrg_curved_center_z_interpolator(data)
    return nd_bspline_interpolation((data.x_axis,), data.center_z; degree = 1, max_derivative_order_eval = 1)
end

function opencrg_curved_heading_interpolator(data)
    return nd_bspline_interpolation((data.x_axis,), data.center_heading; degree = 1, max_derivative_order_eval = 1)
end

function opencrg_curved_height_interpolator(data)
    return nd_bspline_interpolation((data.x_axis, data.v_axis), data.heights; degree = (1, 1), max_derivative_order_eval = 1)
end

function opencrg_curved_surface_interpolator(data)
    center_z = opencrg_curved_center_z_interpolator(data)
    center_heading = opencrg_curved_heading_interpolator(data)
    height = opencrg_curved_height_interpolator(data)
    return (x, z) -> begin
        heading = center_heading(x)
        lateral = cos(heading) * (z - center_z(x))
        height(x, lateral)
    end
end

country_road_center_z_interpolator(; kwargs...) = opencrg_curved_center_z_interpolator(country_road_curved_data(; kwargs...))
country_road_heading_interpolator(; kwargs...) = opencrg_curved_heading_interpolator(country_road_curved_data(; kwargs...))
country_road_surface_interpolator(; kwargs...) = opencrg_curved_surface_interpolator(country_road_curved_data(; kwargs...))
