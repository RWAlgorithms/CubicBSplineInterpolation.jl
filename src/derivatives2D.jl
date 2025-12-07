# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

#### 2D each based on query2D()

function _eval_basis_splines2D(
        C, k1_lb, k2_lb,
        spline1_0, spline1_1, spline1_2, spline1_3,
        spline2_0, spline2_1, spline2_2, spline2_3,
    )

    # The 4x4 = 16 terms.
    out00 = C[begin + k1_lb, begin + k2_lb] * spline1_0 * spline2_0
    out10 = C[begin + k1_lb + 1, begin + k2_lb] * spline1_1 * spline2_0
    out20 = C[begin + k1_lb + 2, begin + k2_lb] * spline1_2 * spline2_0
    out30 = C[begin + k1_lb + 3, begin + k2_lb] * spline1_3 * spline2_0

    out01 = C[begin + k1_lb, begin + k2_lb + 1] * spline1_0 * spline2_1
    out11 = C[begin + k1_lb + 1, begin + k2_lb + 1] * spline1_1 * spline2_1
    out21 = C[begin + k1_lb + 2, begin + k2_lb + 1] * spline1_2 * spline2_1
    out31 = C[begin + k1_lb + 3, begin + k2_lb + 1] * spline1_3 * spline2_1

    out02 = C[begin + k1_lb, begin + k2_lb + 2] * spline1_0 * spline2_2
    out12 = C[begin + k1_lb + 1, begin + k2_lb + 2] * spline1_1 * spline2_2
    out22 = C[begin + k1_lb + 2, begin + k2_lb + 2] * spline1_2 * spline2_2
    out32 = C[begin + k1_lb + 3, begin + k2_lb + 2] * spline1_3 * spline2_2

    out03 = C[begin + k1_lb, begin + k2_lb + 3] * spline1_0 * spline2_3
    out13 = C[begin + k1_lb + 1, begin + k2_lb + 3] * spline1_1 * spline2_3
    out23 = C[begin + k1_lb + 2, begin + k2_lb + 3] * spline1_2 * spline2_3
    out33 = C[begin + k1_lb + 3, begin + k2_lb + 3] * spline1_3 * spline2_3

    return out00 + out10 + out20 + out30 +
        out01 + out11 + out21 + out31 +
        out02 + out12 + out22 + out32 +
        out03 + out13 + out23 + out33
end

#### first derivatives

function query2D_derivative1(x1_in::T, x2_in::T, itp::Interpolator2D{T}) where {T <: AbstractFloat}
    C, A1, A2 = itp.coeffs, itp.x1_query_cache, itp.x2_query_cache

    # # Transform clamp the input coordinates.
    x1 = to_std_interval(x1_in, A1.a, A1.d_div_bma)
    x2 = to_std_interval(x2_in, A2.a, A2.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)

    x1_ub = T(size(C, 1) - 1 - 2)
    if x1 < x_lb
        x1 = x_lb
    elseif x1 > x1_ub
        x1 = x1_ub
    end

    x2_ub = T(size(C, 2) - 1 - 2)
    if x2 < x_lb
        x2 = x_lb
    elseif x2 > x2_ub
        x2 = x2_ub
    end

    # # Compute the 4x4 = 16 terms.
    k1_lb = trunc(Int, x1) - 1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.
    k2_lb = trunc(Int, x2) - 1

    # wrt dx1, x2 is held constant at itp evaluation.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_dB3(x1 - k1_lb)
    spline1_1 = eval_dB3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_dB3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_dB3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d_itp_dx1 = _eval_basis_splines2D(C, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # wrt dx2, x1 is held constant at itp evaluation.
    spline2_0 = eval_dB3(x2 - k2_lb)
    spline2_1 = eval_dB3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_dB3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_dB3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d_itp_dx2 = _eval_basis_splines2D(C, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # apply chain rule.
    d_x1 = A1.d_div_bma
    d_x2 = A2.d_div_bma
    return d_itp_dx1 * d_x1, d_itp_dx2 * d_x2
end

function query2D_derivative1(x1_in::T, x2_in::T, itp::Interpolator2DComplex{T}) where {T <: AbstractFloat}

    Cr, Ci, A1, A2 = itp.real_coeffs, itp.imag_coeffs, itp.x1_query_cache, itp.x2_query_cache

    ### copied from the other query2D
    # Not re-using code from the other query2D  since I want the entire algorithm in a single scope so that the compiler might be able to optimize things better in the same scope.

    # # Transform clamp the input coordinates.
    x1 = to_std_interval(x1_in, A1.a, A1.d_div_bma)
    x2 = to_std_interval(x2_in, A2.a, A2.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)

    x1_ub = T(size(Cr, 1) - 1 - 2)
    if x1 < x_lb
        x1 = x_lb
    elseif x1 > x1_ub
        x1 = x1_ub
    end

    x2_ub = T(size(Cr, 2) - 1 - 2)
    if x2 < x_lb
        x2 = x_lb
    elseif x2 > x2_ub
        x2 = x2_ub
    end

    # # Compute the 4x4 = 16 terms.
    k1_lb = trunc(Int, x1) - 1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.
    k2_lb = trunc(Int, x2) - 1

    # wrt dx1, x2 is held constant at itp evaluation.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_dB3(x1 - k1_lb)
    spline1_1 = eval_dB3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_dB3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_dB3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d_itp_dx1_real = _eval_basis_splines2D(Cr, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)
    d_itp_dx1_imag = _eval_basis_splines2D(Ci, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # wrt dx2, x1 is held constant at itp evaluation.
    spline2_0 = eval_dB3(x2 - k2_lb)
    spline2_1 = eval_dB3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_dB3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_dB3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d_itp_dx2_real = _eval_basis_splines2D(Cr, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)
    d_itp_dx2_imag = _eval_basis_splines2D(Ci, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # apply chain rule.
    d_x1 = A1.d_div_bma
    d_x2 = A2.d_div_bma

    # derivatives wrt x1 has 2 terms, derivatives wrt x2 has 2 terms.
    return d_itp_dx1_real * d_x1, d_itp_dx2_real * d_x2,
        d_itp_dx1_imag * d_x1, d_itp_dx2_imag * d_x2
end

#### second derivatives

function query2D_derivative2(x1_in::T, x2_in::T, itp::Interpolator2D{T}) where {T <: AbstractFloat}
    C, A1, A2 = itp.coeffs, itp.x1_query_cache, itp.x2_query_cache

    # # Transform clamp the input coordinates.
    x1 = to_std_interval(x1_in, A1.a, A1.d_div_bma)
    x2 = to_std_interval(x2_in, A2.a, A2.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)

    x1_ub = T(size(C, 1) - 1 - 2)
    if x1 < x_lb
        x1 = x_lb
    elseif x1 > x1_ub
        x1 = x1_ub
    end

    x2_ub = T(size(C, 2) - 1 - 2)
    if x2 < x_lb
        x2 = x_lb
    elseif x2 > x2_ub
        x2 = x2_ub
    end

    # # Compute the 4x4 = 16 terms.
    k1_lb = trunc(Int, x1) - 1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.
    k2_lb = trunc(Int, x2) - 1

    # wrt dx1, order 2.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_d2B3(x1 - k1_lb)
    spline1_1 = eval_d2B3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_d2B3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_d2B3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx1dx1 = _eval_basis_splines2D(C, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # wrt dx1 dx2, order 1 each.
    spline2_0 = eval_dB3(x2 - k2_lb)
    spline2_1 = eval_dB3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_dB3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_dB3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_dB3(x1 - k1_lb)
    spline1_1 = eval_dB3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_dB3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_dB3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx1dx2 = _eval_basis_splines2D(C, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # wrt dx2, order 2.
    spline2_0 = eval_d2B3(x2 - k2_lb)
    spline2_1 = eval_d2B3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_d2B3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_d2B3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx2dx2 = _eval_basis_splines2D(C, k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3)

    # apply chain rule.
    d_x1 = A1.d_div_bma
    d_x2 = A2.d_div_bma
    return d2_itp_dx1dx1 * d_x1^2, d2_itp_dx1dx2 * d_x1 * d_x2, d2_itp_dx2dx2 * d_x2^2
end

function query2D_derivative2(x1_in::T, x2_in::T, itp::Interpolator2DComplex{T}) where {T <: AbstractFloat}

    Cr, Ci, A1, A2 = itp.real_coeffs, itp.imag_coeffs, itp.x1_query_cache, itp.x2_query_cache

    ### copied from the other query2D
    # Not re-using code from the other query2D  since I want the entire algorithm in a single scope so that the compiler might be able to optimize things better in the same scope.

    # # Transform clamp the input coordinates.
    x1 = to_std_interval(x1_in, A1.a, A1.d_div_bma)
    x2 = to_std_interval(x2_in, A2.a, A2.d_div_bma)

    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)

    x1_ub = T(size(Cr, 1) - 1 - 2)
    if x1 < x_lb
        x1 = x_lb
    elseif x1 > x1_ub
        x1 = x1_ub
    end

    x2_ub = T(size(Cr, 2) - 1 - 2)
    if x2 < x_lb
        x2 = x_lb
    elseif x2 > x2_ub
        x2 = x2_ub
    end

    # # Compute the 4x4 = 16 terms.
    k1_lb = trunc(Int, x1) - 1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.
    k2_lb = trunc(Int, x2) - 1

    # wrt dx1, order 2.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_d2B3(x1 - k1_lb)
    spline1_1 = eval_d2B3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_d2B3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_d2B3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx1dx1_real = _eval_basis_splines2D(
        Cr,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )
    d2_itp_dx1dx1_imag = _eval_basis_splines2D(
        Ci,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )

    # wrt dx1 dx2, order 1 each.
    spline2_0 = eval_dB3(x2 - k2_lb)
    spline2_1 = eval_dB3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_dB3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_dB3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_dB3(x1 - k1_lb)
    spline1_1 = eval_dB3(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_dB3(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_dB3(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx1dx2_real = _eval_basis_splines2D(
        Cr,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )
    d2_itp_dx1dx2_imag = _eval_basis_splines2D(
        Ci,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )

    # wrt dx2, order 2.
    spline2_0 = eval_d2B3(x2 - k2_lb)
    spline2_1 = eval_d2B3(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_d2B3(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_d2B3(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb - 3) # x1 - (k1_lb + 3)

    d2_itp_dx2dx2_real = _eval_basis_splines2D(
        Cr,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )
    d2_itp_dx2dx2_imag = _eval_basis_splines2D(
        Ci,
        k1_lb, k2_lb, spline1_0, spline1_1, spline1_2, spline1_3, spline2_0, spline2_1, spline2_2, spline2_3,
    )

    # apply chain rule.
    d_x1 = A1.d_div_bma
    d_x2 = A2.d_div_bma

    # real part has 3 terms, imag part has 3 terms.
    return d2_itp_dx1dx1_real * d_x1^2, d2_itp_dx1dx2_real * d_x1 * d_x2, d2_itp_dx2dx2_real * d_x2^2,
        d2_itp_dx1dx1_imag * d_x1^2, d2_itp_dx1dx2_imag * d_x1 * d_x2, d2_itp_dx2dx2_imag * d_x2^2
end
