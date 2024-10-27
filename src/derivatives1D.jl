# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>



# TODO unroll and call specialized versions to reduce branch prediction.
function eval_dB3(x::T) where T <: AbstractFloat
    return eval_B2(x + half(T)) - eval_B2(x - half(T))
end

function eval_d2B3(x::T) where T <: AbstractFloat
    return eval_dB2(x + half(T)) - eval_dB2(x - half(T))
end

function eval_dB2(x::T) where T <: AbstractFloat
    return eval_B1(x + half(T)) - eval_B1(x - half(T))
end

# TODO specialize according to region
function eval_B2(x::T) where T <: AbstractFloat
    
    abs_x = abs(x)
    if abs_x < half(T)
        return T(3)/T(4) -x^2
    elseif half(T) <= abs_x < T(3)/T(2)
        return ((3 - 2*abs_x)^2)/8
    end

    return zero(T)
end

# TODO specialize according to region
function eval_B1(x::T) where T <: AbstractFloat
    
    abs_x = abs(x)
    if zero(T) < abs_x < one(T)
        return one(T) - abs_x
    end

    return zero(T)
end

#### 1D

# based on query1D
function query1D_derivative1(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    x_lb = T(2)
    
    x_ub = T(size(c,1)-1-2)
    if x < x_lb
        x = x_lb
    elseif x > x_ub
        x = x_ub
    end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x)-1
    
    k = k_lb
    #out1 = c[begin + k]*eval_cubic_spline(x - k) # non-specialized.
    #out1 = c[begin + k]*eval_dB3_in12(x - k)
    out1 = c[begin + k]*eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in12(x - k))

    k = k_lb + 1
    #out2 = c[begin + k]*eval_dB3_in01(x - k)
    out2 = c[begin + k]*eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in01(x - k))

    k = k_lb + 2
    #out3 = c[begin + k]*eval_dB3_in01(x - k)
    out3 = c[begin + k]*eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in01(x - k))
    
    k = k_lb + 3
    #out4 = c[begin + k]*eval_dB3_in12(x - k)
    out4 = c[begin + k]*eval_dB3(x - k)
    #@assert isapprox(eval_dB3(x - k), eval_dB3_in12(x - k))

    wrt_x = out1 + out2 + out3 + out4
    d_x_wrt_x_in = A.d_div_bma
    return wrt_x*d_x_wrt_x_in
end

# based on query1D
function query1D_derivative1(x_in::T, itp::Interpolator1DComplex{T}) where T <: AbstractFloat
    
    cr, ci, A = itp.real_coeffs, itp.imag_coeffs, itp.query_cache
    
    ### copied from the other query1D
    # Not re-using code from the other query2D  since I want the entire algorithm in a single scope so that the compiler might be able to optimize things better in the same scope.

    # # Transform clamp the input coordinates.
    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)
    
    x_ub = T(size(cr,1)-1-2)
    if x < x_lb
        x = x_lb
    elseif x > x_ub
        x = x_ub
    end

    # # Compute the 4x4 = 16 terms.
    k_lb = trunc(Int, x)-1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.

    # Compute shared terms spline1_2 means direction x, shift by 2 from the lower bound of the sampling lattice.
    # the point x is near the center of the 4x4 window on the sampling lattice, so shift 2 and 3 use the 01 version of the spline, similarly for shifts 0 and 4 for the 12 version.
    spline_0 = eval_dB3(x - k_lb)
    spline_1 = eval_dB3(x - k_lb - 1) # x - (k_lb + 1)
    spline_2 = eval_dB3(x - k_lb - 2) # x - (k_lb + 2)
    spline_3 = eval_dB3(x - k_lb - 3) # x - (k_lb + 3)
    #### end copy from the other query1D.

    # real.
    out_r_0 = cr[begin + k_lb]*spline_0
    out_r_1 = cr[begin + k_lb + 1]*spline_1
    out_r_2 = cr[begin + k_lb + 2]*spline_2
    out_r_3 = cr[begin + k_lb + 3]*spline_3

    out_real = out_r_0 + out_r_1 + out_r_2 + out_r_3
    
    # imaginary.
    out_i_0 = ci[begin + k_lb]*spline_0
    out_i_1 = ci[begin + k_lb + 1]*spline_1
    out_i_2 = ci[begin + k_lb + 2]*spline_2
    out_i_3 = ci[begin + k_lb + 3]*spline_3

    out_imag = out_i_0 + out_i_1 + out_i_2 + out_i_3
    
    d_x_wrt_x_in = A.d_div_bma
    return out_real*d_x_wrt_x_in, out_imag*d_x_wrt_x_in
end

# # Second derivative

# based on query1D
function query1D_derivative2(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    x_lb = T(2)
    
    x_ub = T(size(c,1)-1-2)
    if x < x_lb
        x = x_lb
    elseif x > x_ub
        x = x_ub
    end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x)-1

    k = k_lb
    out1 = c[begin + k]*eval_d2B3(x - k)
    #@assert 1 <= abs(x - k) <= 2

    k = k_lb + 1
    out2 = c[begin + k]*eval_d2B3(x - k)
    #@assert 0 <= abs(x - k) <= 1

    k = k_lb + 2
    out3 = c[begin + k]*eval_d2B3(x - k)
    #@assert 0 <= abs(x - k) <= 1
    
    k = k_lb + 3
    out4 = c[begin + k]*eval_d2B3(x - k)
    #@assert 1 <= abs(x - k) <= 2

    
    wrt_x = out1 + out2 + out3 + out4
    d_x_wrt_x_in = A.d_div_bma
    return wrt_x*d_x_wrt_x_in^2
end

# based on query1D
function query1D_derivative2(x_in::T, itp::Interpolator1DComplex{T}) where T <: AbstractFloat
    
    cr, ci, A = itp.real_coeffs, itp.imag_coeffs, itp.query_cache
    
    ### copied from the other query1D
    # Not re-using code from the other query2D  since I want the entire algorithm in a single scope so that the compiler might be able to optimize things better in the same scope.

    # # Transform clamp the input coordinates.
    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.

    # manual if statement is faster than clamp for some reason.
    x_lb = T(2)
    
    x_ub = T(size(cr,1)-1-2)
    if x < x_lb
        x = x_lb
    elseif x > x_ub
        x = x_ub
    end

    # # Compute the 4x4 = 16 terms.
    k_lb = trunc(Int, x)-1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.

    # Compute shared terms spline1_2 means direction x, shift by 2 from the lower bound of the sampling lattice.
    # the point x is near the center of the 4x4 window on the sampling lattice, so shift 2 and 3 use the 01 version of the spline, similarly for shifts 0 and 4 for the 12 version.
    spline_0 = eval_d2B3(x - k_lb)
    spline_1 = eval_d2B3(x - k_lb - 1) # x - (k_lb + 1)
    spline_2 = eval_d2B3(x - k_lb - 2) # x - (k_lb + 2)
    spline_3 = eval_d2B3(x - k_lb - 3) # x - (k_lb + 3)
    #### end copy from the other query1D.

    # real.
    out_r_0 = cr[begin + k_lb]*spline_0
    out_r_1 = cr[begin + k_lb + 1]*spline_1
    out_r_2 = cr[begin + k_lb + 2]*spline_2
    out_r_3 = cr[begin + k_lb + 3]*spline_3

    out_real = out_r_0 + out_r_1 + out_r_2 + out_r_3
    
    # imaginary.
    out_i_0 = ci[begin + k_lb]*spline_0
    out_i_1 = ci[begin + k_lb + 1]*spline_1
    out_i_2 = ci[begin + k_lb + 2]*spline_2
    out_i_3 = ci[begin + k_lb + 3]*spline_3

    out_imag = out_i_0 + out_i_1 + out_i_2 + out_i_3
    
    d_x_wrt_x_in_squared = A.d_div_bma^2
    return out_real*d_x_wrt_x_in_squared, out_imag*d_x_wrt_x_in_squared
end
