# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>


# # Complex-valued 1D

struct Interpolator1DComplex{T <: AbstractFloat} <: AbstractInterpolator1D
    
    real_coeffs::Memory{T}
    imag_coeffs::Memory{T}
    query_cache::IntervalConversion{T}

    x_start::T
    x_fin::T

    # # the boundary 2 samples won't won't the provided samples.
    # function Interpolator1DComplex(s_real::Union{Memory{T},Vector{T}}, s_imag::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
    #     size(s_real) == size(s_imag) || error("Size mismatch between the real and imaginary sample matrices.")
        
    #     A = IntervalConversion(x_start, x_fin, length(s_real))
        
    #     c_real = Memory{T}(undef, length(s_real))
    #     get_coeffs!(c_real, s_real, ϵ)

    #     c_imag = Memory{T}(undef, length(s_imag))
    #     get_coeffs!(c_imag, s_imag, ϵ)
        
    #     # sanity check. The following should pass since S_real and S_imag have the same size.
    #     # Keep it here to explicitly define invariants on real_coeffs and image_coeffs.
    #     size(c_real) == size(c_imag) || error("Size mismatch for coefficients. Please file bug report.")
    #     return new{T}(c_real, c_imag, A, x_start, x_fin)
    # end

    # Mutates buf, option is for dispatch.
    function Interpolator1DComplex(
        padding_option::PaddingOption,
        extrapolation_option::ExtrapolationOption,
        buf::FitBuffer1D,
        s_real::Union{Memory{T},Vector{T}},
        s_imag::Union{Memory{T},Vector{T}},
        x_start::T,
        x_fin::T;
        ϵ::T = eps(T)*2,
        ) where T <: AbstractFloat

        size(s_real) == size(s_imag) || error("Size mismatch between the real and imaginary sample matrices.")

        x_start < x_fin || error("x_start must be strictly smaller than x_fin.")

        Np = buf.N_padding

        A, tmp_r = create_query_cache(x_start, x_fin, length(s_real), Np)
        
        c_real = Memory{T}(undef, get_num_coeffs(buf))
        get_coeffs!(padding_option, extrapolation_option, c_real, buf, s_real, tmp_r, ϵ)

        c_imag = Memory{T}(undef, get_num_coeffs(buf))
        get_coeffs!(padding_option, extrapolation_option, c_imag, buf, s_imag, tmp_r, ϵ)

        size(c_real) == size(c_imag) || error("Size mismatch for coefficients. Please file bug report.")
        return new{T}(c_real, c_imag, A, x_start, x_fin)
    end

    # convinence constructor. Mutates buf
    function Interpolator1DComplex(buf::FitBuffer1D, s_real::Union{Memory{T},Vector{T}}, s_imag::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
        return Interpolator1DComplex(LinearPadding(), ConstantExtrapolation(), buf, s_real, s_imag, x_start, x_fin; ϵ = ϵ)
    end    
end

# Mutates buf, option is for dispatch. itp mutates, is output.
function update_itp!(
    padding_option::PaddingOption,
    extrapolation_option::ExtrapolationOption,
    itp::Interpolator1DComplex,
    buf::FitBuffer1D,
    s_real::Union{Memory{T},Vector{T}},
    s_imag::Union{Memory{T},Vector{T}};
    ϵ::T = eps(T)*2,
    ) where T <: AbstractFloat

    length(s_real) == length(s_imag) == get_data_length(buf) || error("Length mismatch.")

    x_start, x_fin = get_itp_interval(itp)
    get_coeffs!(padding_option, extrapolation_option, itp.real_coeffs, buf, s_real, LinRange(x_start, x_fin, length(s_real)), ϵ)
    get_coeffs!(padding_option, extrapolation_option, itp.imag_coeffs, buf, s_imag, LinRange(x_start, x_fin, length(s_imag)), ϵ)
    return nothing
end

# convenince
function update_itp!(itp::Interpolator1DComplex, buf::FitBuffer1D, s_real::Union{Memory{T},Vector{T}}, s_imag::Union{Memory{T},Vector{T}}; ϵ::T = eps(T)*2) where T <: AbstractFloat
    return update_itp!(LinearPadding(), ConstantExtrapolation(), itp, buf, s_real, s_imag; ϵ = ϵ)
end


function query1D(x_in::T, itp::Interpolator1DComplex{T}) where T <: AbstractFloat
    
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
    spline_0 = eval_cubic_spline_in12(x - k_lb)
    spline_1 = eval_cubic_spline_in01(x - k_lb - 1) # x - (k_lb + 1)
    spline_2 = eval_cubic_spline_in01(x - k_lb - 2) # x - (k_lb + 2)
    spline_3 = eval_cubic_spline_in12(x - k_lb - 3) # x - (k_lb + 3)
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
    
    return Complex(out_real, out_imag)
end


#### 2D

struct Interpolator2DComplex{T <: AbstractFloat} <: AbstractInterpolator2D
    real_coeffs::Matrix{T}
    imag_coeffs::Matrix{T}
    x1_query_cache::IntervalConversion{T}
    x2_query_cache::IntervalConversion{T}

    x1_start::T
    x1_fin::T
    x2_start::T
    x2_fin::T

    # # no padding. The bounary two samples will not patch the provided samples.
    # function Interpolator2DComplex(S_real::Matrix{T}, S_imag::Matrix{T}, x1_start::T, x1_fin::T, x2_start::T, x2_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
    #     size(S_real) == size(S_imag) || error("Size mismatch between the real and imaginary sample matrices.")
        
    #     A1 = IntervalConversion(x1_start, x1_fin, size(S_real,1))
    #     A2 = IntervalConversion(x2_start, x2_fin, size(S_real,2))

    #     X = zeros(T, size(S_real))
    #     buf1 = Memory{T}(undef, size(S_real,1))
    #     buf2 = Memory{T}(undef, size(S_real,2))

    #     Yr = zeros(T, size(S_real))
    #     get_coeffs!(Yr, X, buf1, buf2, S_real, ϵ)

    #     Yi = zeros(T, size(S_imag))
    #     get_coeffs!(Yi, X, buf1, buf2, S_real, ϵ)

    #     return new{T}(Yr, Yi, A1, A2, )
    # end

    # has padding.
    # mutates buf.
    function Interpolator2DComplex(
        padding_option::PaddingOption,
        extrapolation_option::ExtrapolationOption,
        buf::FitBuffer2D,
        S_real::Matrix{T},
        S_imag::Matrix{T},
        x1_start::T,
        x1_fin::T,
        x2_start::T,
        x2_fin::T;
        ϵ::T = eps(T)*2,
        ) where T <: AbstractFloat

        size(S_real) == size(S_imag) || error("Size mismatch between the real and imaginary sample matrices.")

        A1, xs1 = create_query_cache(x1_start, x1_fin, size(S_real,1), buf.buf_x1.N_padding)
        A2, xs2 = create_query_cache(x2_start, x2_fin, size(S_real,2), buf.buf_x2.N_padding)

        c_real = zeros(T, get_coeffs_size(buf)) # allocates.
        get_coeffs!(padding_option, extrapolation_option, c_real, buf, S_real, xs1, xs2, ϵ)

        c_imag = zeros(T, get_coeffs_size(buf)) # allocates.
        get_coeffs!(padding_option, extrapolation_option, c_imag, buf, S_imag, xs1, xs2, ϵ)
    
        return new{T}(c_real, c_imag, A1, A2, x1_start, x1_fin, x2_start, x2_fin)
    end

    # has padding. Default to LinearPadding()
    function Interpolator2DComplex(buf::FitBuffer2D, S_real::Matrix{T}, S_imag::Matrix{T}, x1_start::T, x1_fin::T, x2_start::T, x2_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
        return Interpolator2DComplex(LinearPadding(), ConstantExtrapolation(), buf,  S_real, S_imag, x1_start, x1_fin, x2_start, x2_fin; ϵ = ϵ)
    end
end

# Mutates buf, option is for dispatch. itp mutates, is output.
function update_itp!(
    padding_option::PaddingOption,
    extrapolation_option::ExtrapolationOption,
    itp::Interpolator2DComplex,
    buf::FitBuffer2D,
    S_real::Matrix{T},
    S_imag::Matrix{T};
    ϵ::T = eps(T)*2,
    ) where T <: AbstractFloat
    
    size(S_real) == size(S_imag) || error("Size mismatch between the real and imaginary sample matrices.")
    size(itp.real_coeffs) == size(itp.imag_coeffs) == get_coeffs_size(buf) || error("Size mismatch.")
    size(S_real,2) == size(buf.S1, 2) || error("Size mismatch.")

    x1_start, x1_fin, x2_start, x2_fin = get_itp_interval(itp)
    xs1 = LinRange(x1_start, x1_fin, length(S_real))
    xs2 = LinRange(x2_start, x2_fin, length(S_real))

    get_coeffs!(padding_option, extrapolation_option, itp.real_coeffs, buf, S_real, xs1, xs2, ϵ)
    get_coeffs!(padding_option, extrapolation_option, itp.imag_coeffs, buf, S_imag, xs1, xs2, ϵ)
    return nothing
end

# convenince
function update_itp!(itp::Interpolator2DComplex, buf::FitBuffer2D, S_real::Matrix{T}, S_imag::Matrix{T}; ϵ::T = eps(T)*2) where T <: AbstractFloat
    return update_itp!(LinearPadding(), ConstantExtrapolation(), itp, buf, S_real, S_imag; ϵ = ϵ)
end

function query2D(x1_in::T, x2_in::T, itp::Interpolator2DComplex{T}) where T <: AbstractFloat
    
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
    
    x1_ub = T(size(Cr,1)-1-2)
    if x1 < x_lb
        x1 = x_lb
    elseif x1 > x1_ub
        x1 = x1_ub
    end
    
    x2_ub = T(size(Cr,2)-1-2)
    if x2 < x_lb
        x2 = x_lb
    elseif x2 > x2_ub
        x2 = x2_ub
    end

    # # Compute the 4x4 = 16 terms.
    k1_lb = trunc(Int, x1)-1 # trunc(Int, .) is a bit faster than ceil(Int, .) in v1.11-rc1, but still takes a long time.
    k2_lb = trunc(Int, x2)-1

    # Compute shared terms spline1_2 means direction x1, shift by 2 from the lower bound of the sampling lattice.
    # the point (x1, x2) is near the center of the 4x4 window on the sampling lattice, so shift 2 and 3 use the 01 version of the spline, similarly for shifts 0 and 4 for the 12 version.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb - 1) # x2 - (k2_lb + 1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb - 2) # x2 - (k2_lb + 2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb - 3) # x2 - (k2_lb + 3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb - 1) # x1 - (k1_lb + 1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb - 2) # x1 - (k1_lb + 2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb - 3) # x1 - (k1_lb + 3)
    #### end copy from the other query2D.

    # cache the tensor product basis
    sp_00 = spline1_0*spline2_0
    sp_10 = spline1_1*spline2_0
    sp_20 = spline1_2*spline2_0
    sp_30 = spline1_3*spline2_0

    sp_01 = spline1_0*spline2_1
    sp_11 = spline1_1*spline2_1
    sp_21 = spline1_2*spline2_1
    sp_31 = spline1_3*spline2_1

    sp_02 = spline1_0*spline2_2
    sp_12 = spline1_1*spline2_2
    sp_22 = spline1_2*spline2_2
    sp_32 = spline1_3*spline2_2

    sp_03 = spline1_0*spline2_3
    sp_13 = spline1_1*spline2_3
    sp_23 = spline1_2*spline2_3
    sp_33 = spline1_3*spline2_3

    # The 4x4 = 16 terms.
    out_r_00 = Cr[begin + k1_lb, begin + k2_lb]*sp_00
    out_r_10 = Cr[begin + k1_lb + 1, begin + k2_lb]*sp_10
    out_r_20 = Cr[begin + k1_lb + 2, begin + k2_lb]*sp_20
    out_r_30 = Cr[begin + k1_lb + 3, begin + k2_lb]*sp_30

    out_r_01 = Cr[begin + k1_lb, begin + k2_lb + 1]*sp_01
    out_r_11 = Cr[begin + k1_lb + 1, begin + k2_lb + 1]*sp_11
    out_r_21 = Cr[begin + k1_lb + 2, begin + k2_lb + 1]*sp_21
    out_r_31 = Cr[begin + k1_lb + 3, begin + k2_lb + 1]*sp_31

    out_r_02 = Cr[begin + k1_lb, begin + k2_lb + 2]*sp_02
    out_r_12 = Cr[begin + k1_lb + 1, begin + k2_lb + 2]*sp_12
    out_r_22 = Cr[begin + k1_lb + 2, begin + k2_lb + 2]*sp_22
    out_r_32 = Cr[begin + k1_lb + 3, begin + k2_lb + 2]*sp_32

    out_r_03 = Cr[begin + k1_lb, begin + k2_lb + 3]*sp_03
    out_r_13 = Cr[begin + k1_lb + 1, begin + k2_lb + 3]*sp_13
    out_r_23 = Cr[begin + k1_lb + 2, begin + k2_lb + 3]*sp_23
    out_r_33 = Cr[begin + k1_lb + 3, begin + k2_lb + 3]*sp_33

    out_real = out_r_00 + out_r_10 + out_r_20 + out_r_30 +
        out_r_01 + out_r_11 + out_r_21 + out_r_31 +
        out_r_02 + out_r_12 + out_r_22 + out_r_32 +
        out_r_03 + out_r_13 + out_r_23 + out_r_33
        
    #
    out_i_00 = Ci[begin + k1_lb, begin + k2_lb]*sp_00
    out_i_10 = Ci[begin + k1_lb + 1, begin + k2_lb]*sp_10
    out_i_20 = Ci[begin + k1_lb + 2, begin + k2_lb]*sp_20
    out_i_30 = Ci[begin + k1_lb + 3, begin + k2_lb]*sp_30

    out_i_01 = Ci[begin + k1_lb, begin + k2_lb + 1]*sp_01
    out_i_11 = Ci[begin + k1_lb + 1, begin + k2_lb + 1]*sp_11
    out_i_21 = Ci[begin + k1_lb + 2, begin + k2_lb + 1]*sp_21
    out_i_31 = Ci[begin + k1_lb + 3, begin + k2_lb + 1]*sp_31

    out_i_02 = Ci[begin + k1_lb, begin + k2_lb + 2]*sp_02
    out_i_12 = Ci[begin + k1_lb + 1, begin + k2_lb + 2]*sp_12
    out_i_22 = Ci[begin + k1_lb + 2, begin + k2_lb + 2]*sp_22
    out_i_32 = Ci[begin + k1_lb + 3, begin + k2_lb + 2]*sp_32

    out_i_03 = Ci[begin + k1_lb, begin + k2_lb + 3]*sp_03
    out_i_13 = Ci[begin + k1_lb + 1, begin + k2_lb + 3]*sp_13
    out_i_23 = Ci[begin + k1_lb + 2, begin + k2_lb + 3]*sp_23
    out_i_33 = Ci[begin + k1_lb + 3, begin + k2_lb + 3]*sp_33

    out_imag = out_i_00 + out_i_10 + out_i_20 + out_i_30 +
        out_i_01 + out_i_11 + out_i_21 + out_i_31 +
        out_i_02 + out_i_12 + out_i_22 + out_i_32 +
        out_i_03 + out_i_13 + out_i_23 + out_i_33
    
    return Complex(out_real, out_imag)
end
