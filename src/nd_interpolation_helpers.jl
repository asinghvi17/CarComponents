using DataInterpolationsND

"""
    bspline_knots_for_control_axis(axis, degree)

Construct an open/clamped knot vector for `DataInterpolationsND.BSplineInterpolationDimension`
from a vector of control-point coordinates. For `n` control points and spline
`degree = p`, DataInterpolationsND expects `n = length(knots_all) - p - 1` basis
functions. With the package's default endpoint multiplicities, this means the
unique knot vector should have length `n - p + 1`.
"""
function bspline_knots_for_control_axis(axis::AbstractVector{<:Real}, degree::Integer)
    n = length(axis)
    @assert degree >= 0 "B-spline degree must be non-negative."
    @assert n >= degree + 1 "Need at least degree + 1 control points. Got n=$(n), degree=$(degree)."
    @assert all(>(0), diff(axis)) "Control-point coordinate axis must be strictly increasing."

    n_unique_knots = n - degree + 1
    return collect(range(first(axis), last(axis); length = n_unique_knots))
end

function _degree_tuple(degree::Integer, ::Val{N}) where {N}
    return ntuple(_ -> Int(degree), Val(N))
end

function _degree_tuple(degree, ::Val{N}) where {N}
    @assert length(degree) == N "Degree tuple length must match number of interpolation dimensions."
    return ntuple(i -> Int(degree[i]), Val(N))
end

"""
    nd_bspline_interpolation(axes, values; degree = 2, max_derivative_order_eval = 1)

Build a scalar `DataInterpolationsND.NDInterpolation` using B-spline dimensions.
`axes` is an N-tuple of strictly increasing control-point coordinate vectors.
`values` is the N-dimensional control lattice. Its first N dimensions must match
`length.(axes)`.
"""
function nd_bspline_interpolation(
        axes::NTuple{N, <:AbstractVector{<:Real}},
        values::AbstractArray;
        degree = 2,
        max_derivative_order_eval::Integer = 1,
    ) where {N}
    degrees = _degree_tuple(degree, Val(N))
    expected_size = ntuple(i -> length(axes[i]), Val(N))
    @assert size(values)[1:N] == expected_size "The first N dimensions of values must match length.(axes). Expected $(expected_size), got $(size(values)[1:N])."

    interp_dims = ntuple(Val(N)) do i
        knots = bspline_knots_for_control_axis(collect(axes[i]), degrees[i])
        DataInterpolationsND.BSplineInterpolationDimension(
            knots,
            degrees[i];
            max_derivative_order_eval = max_derivative_order_eval,
        )
    end

    return DataInterpolationsND.NDInterpolation(values, interp_dims)
end

function _validate_rectilinear_surface_matrices(X::AbstractMatrix, Y::AbstractMatrix, Z::AbstractMatrix)
    @assert size(X) == size(Y) == size(Z) "X, Y, and Z matrices must have the same size."
    x_axis = collect(X[:, 1])
    y_axis = collect(Y[1, :])
    @assert all(>(0), diff(x_axis)) "X[:, 1] must be strictly increasing."
    @assert all(>(0), diff(y_axis)) "Y[1, :] must be strictly increasing."

    for j in axes(X, 2)
        @assert all(isapprox.(X[:, j], x_axis; rtol = 1e-10, atol = 1e-12)) "X must be rectilinear: every column must match X[:, 1]."
    end
    for i in axes(Y, 1)
        @assert all(isapprox.(Y[i, :], y_axis; rtol = 1e-10, atol = 1e-12)) "Y must be rectilinear: every row must match Y[1, :]."
    end

    return x_axis, y_axis
end

"""
    nd_bspline_interpolation_2d_from_matrices(X, Y, Z, degree_x, degree_y, max_derivative_order_eval)

Build a 2D scalar B-spline interpolator from rectilinear `X`, `Y`, and `Z`
matrices. `X` and `Y` define the control-point coordinate lattice, and `Z` is the
scalar control lattice.
"""
function nd_bspline_interpolation_2d_from_matrices(
        X::AbstractMatrix{<:Real},
        Y::AbstractMatrix{<:Real},
        Z::AbstractMatrix{<:Real},
        degree_x::Integer = 2,
        degree_y::Integer = 2,
        max_derivative_order_eval::Integer = 1,
    )
    x_axis, y_axis = _validate_rectilinear_surface_matrices(X, Y, Z)
    return nd_bspline_interpolation(
        (x_axis, y_axis),
        Z;
        degree = (degree_x, degree_y),
        max_derivative_order_eval = max_derivative_order_eval,
    )
end

function fake_bspline_surface_X(; nx::Integer = 6, ny::Integer = 7)
    xs = collect(range(-1.0, 1.0; length = nx))
    ys = collect(range(-1.5, 1.5; length = ny))
    return [x for x in xs, y in ys]
end

function fake_bspline_surface_Y(; nx::Integer = 6, ny::Integer = 7)
    xs = collect(range(-1.0, 1.0; length = nx))
    ys = collect(range(-1.5, 1.5; length = ny))
    return [y for x in xs, y in ys]
end

function fake_bspline_surface_Z(; nx::Integer = 6, ny::Integer = 7)
    X = fake_bspline_surface_X(; nx, ny)
    Y = fake_bspline_surface_Y(; nx, ny)
    return @. 0.4 + 0.2 * X - 0.1 * Y + 0.05 * X * Y + 0.03 * X^2
end

function fake_bspline_surface_interpolator()
    return nd_bspline_interpolation_2d_from_matrices(
        fake_bspline_surface_X(),
        fake_bspline_surface_Y(),
        fake_bspline_surface_Z(),
        2,
        2,
        1,
    )
end

fake_bspline_surface_value(x, y) = fake_bspline_surface_interpolator()(x, y)
