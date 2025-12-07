# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# This file contain methods that utilize `Interpolator1D{T}` as a minimalist spline interpolator, where the coefficients are externally determined. Extrapolation returns `T(NaN)`. The constructor to use is: `Interpolator1D(lb::T, ub::T, N::Int)`.

function _get_coeffs_range(lb::Real, ub::Real, N::Integer)
    Δx = (ub - lb) / (N - 1)
    return LinRange(lb - 3 * Δx, ub + 3 * Δx, N) # +/- 3 times instead of 2 to ensure queries at lb and ub don't return NaN.
end

"""
    get_coeffs_range(itp::Interpolator1D)

Returns a `LinRange` type where the `n`-th element is the input space position for the `n`-th coefficient in `itp`. This can be used to initialize the coefficients to approximate some callable scalar-valued 1-D function `f` by:

```
itp.coeffs .= f.(get_coeffs_range(itp))
```
"""
function get_coeffs_range(itp::Interpolator1D)
    return _get_coeffs_range(itp.x_start, itp.x_fin, length(itp.coeffs))
end

"""
    get_query_bounds(itp::Interpolator1D{T})

Returns the lower bound (lb) and upper bound (ub), each having type `T`. The bounds create an interval `[lb,ub]` such that `itp1D`, `itp1D_derivative1`, and `itp1D_parameter_derivatives`, are allowed. Querying outside of this interval is considered undefined behavior, and the caller of  `itp1D`, `itp1D_derivative1`, or `itp1D_parameter_derivatives` should ensure that does not happen.
"""
function get_query_bounds(itp::Interpolator1D)
    return itp.x_start, itp.x_fin
end

"""
    itp1D(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat

Identical to `query1D`, except extrapolation outside the interval `[lb, ub]`, given by `lb, ub = get_query_bounds(itp)` returns `T(NaN)`.
"""
function itp1D(x_in::T, itp::Interpolator1D{T}) where {T <: AbstractFloat}
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    x_lb = T(2)

    x_ub = T(size(c, 1) - 1 - 2)
    if x < x_lb || x > x_ub
        return T(NaN)
    end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x) - 1

    k = k_lb
    #out1 = c[begin + k]*eval_cubic_spline(x - k) # non-specialized.
    out1 = c[begin + k] * eval_cubic_spline_in12(x - k)

    k = k_lb + 1
    out2 = c[begin + k] * eval_cubic_spline_in01(x - k)

    k = k_lb + 2
    out3 = c[begin + k] * eval_cubic_spline_in01(x - k)

    k = k_lb + 3
    out4 = c[begin + k] * eval_cubic_spline_in12(x - k)

    return out1 + out2 + out3 + out4
end

"""
    itp1D_derivative1(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat

Identical to `query1D_derivative1`, except extrapolation outside the interval `[lb, ub]`, given by `lb, ub = get_query_bounds(itp)` returns `T(NaN)`.
"""
function itp1D_derivative1(x_in::T, itp::Interpolator1D{T}) where {T <: AbstractFloat}
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    x_lb = T(2)

    x_ub = T(size(c, 1) - 1 - 2)
    if x < x_lb || x > x_ub
        return T(NaN)
    end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x) - 1

    k = k_lb
    #out1 = c[begin + k]*eval_cubic_spline(x - k) # non-specialized.
    #out1 = c[begin + k]*eval_dB3_in12(x - k)
    out1 = c[begin + k] * eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in12(x - k))

    k = k_lb + 1
    #out2 = c[begin + k]*eval_dB3_in01(x - k)
    out2 = c[begin + k] * eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in01(x - k))

    k = k_lb + 2
    #out3 = c[begin + k]*eval_dB3_in01(x - k)
    out3 = c[begin + k] * eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in01(x - k))

    k = k_lb + 3
    #out4 = c[begin + k]*eval_dB3_in12(x - k)
    out4 = c[begin + k] * eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in12(x - k))

    wrt_x = out1 + out2 + out3 + out4
    d_x_wrt_x_in = A.d_div_bma
    return wrt_x * d_x_wrt_x_in
end


"""
    itp1D_parameter_derivatives(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat

Identical to `query1D_derivative1`, except extrapolation outside the interval `[lb, ub]`, given by `lb, ub = get_query_bounds(itp)` return `T(NaN)` for evaluation and derivative values, and `0` for index values.

Return positions: interpolation evaluation, parameter index 1, parameter index 2, parameter index 3, parameter index 4, parameter derivative 1, parameter derivative 2, parameter derivative 3, parameter derivative 4.
"""
function itp1D_parameter_derivatives(x_in::T, itp::Interpolator1D{T}) where {T <: AbstractFloat}
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    x_lb = T(2)

    x_ub = T(size(c, 1) - 1 - 2)
    if x < x_lb || x > x_ub
        return T(NaN), 0, 0, 0, 0, T(NaN), T(NaN), T(NaN), T(NaN)
    end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x) - 1

    k1 = k_lb
    #out1 = c[begin + k]*eval_cubic_spline(x - k) # non-specialized.
    d1 = eval_cubic_spline_in12(x - k1)
    out1 = c[begin + k1] * d1

    k2 = k_lb + 1
    d2 = eval_cubic_spline_in01(x - k2)
    out2 = c[begin + k2] * d2

    k3 = k_lb + 2
    d3 = eval_cubic_spline_in01(x - k3)
    out3 = c[begin + k3] * d3

    k4 = k_lb + 3
    d4 = eval_cubic_spline_in12(x - k4)
    out4 = c[begin + k4] * d4

    return out1 + out2 + out3 + out4, k1 + 1, k2 + 1, k3 + 1, k4 + 1, d1, d2, d3, d4
end
