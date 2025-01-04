# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>



# pre-allocated version
struct FitBuffer2D{T<:AbstractFloat}

    S1::Matrix{T} # intermediate coefficient matris.
    #mat_view1::ST # for filtering in the x1 (row) direction.
    #mat_view2::ST # for x2 direction.

    buf_x1::FitBuffer1D{T}
    buf_x2::FitBuffer1D{T}

    function FitBuffer2D(::Type{T}, sz::Tuple{Int,Int}; N_padding::Tuple{Int,Int}=(10, 10)) where {T<:AbstractFloat}

        Np1, Np2 = N_padding
        Np1 >= 5 || error("All entries of N_padding must be larger or equal to 5. Two for the query system used in this library, three for transition to constant extrapolation.")
        Np2 >= 5 || error("All entries of N_padding must be larger or equal to 5. Two for the query system used in this library, three for transition to constant extrapolation.")

        N1, N2 = sz
        N1 > 2 || error("All entries of sz must be larger than 2.")
        N2 > 2 || error("All entries of sz must be larger than 2.")

        S1 = zeros(T, N1 + 2 * Np1, N2)
        return new{T}(
            S1,
            FitBuffer1D(T, N1; N_padding=Np1),
            FitBuffer1D(T, N2; N_padding=Np2),
        )
    end
end

function get_coeffs_size(s::FitBuffer2D)
    M1 = get_num_coeffs(s.buf_x1)
    M2 = get_num_coeffs(s.buf_x2)
    return (M1, M2)
end

struct Interpolator2D{T<:AbstractFloat} <: AbstractInterpolator2D
    coeffs::Matrix{T}
    x1_query_cache::IntervalConversion{T}
    x2_query_cache::IntervalConversion{T}

    x1_start::T
    x1_fin::T
    x2_start::T
    x2_fin::T

    # # no padding.
    # function Interpolator2D(S::Matrix{T}, x1_start::T, x1_fin::T, x2_start::T, x2_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat

    #     A1 = IntervalConversion(x1_start, x1_fin, size(S,1))
    #     A2 = IntervalConversion(x2_start, x2_fin, size(S,2))

    #     Y = zeros(T, size(S))
    #     X = zeros(T, size(S))
    #     buf1 = Memory{T}(undef, size(S,1))
    #     buf2 = Memory{T}(undef, size(S,2))
    #     _get_coeffs!(Y, X, buf1, buf2, S, ϵ)

    #     return new{T}(Y, A1, A2, x1_start, x1_fin, x2_start, x2_fin)
    # end

    # has padding.
    # mutates buf.
    function Interpolator2D(
        padding_option::PaddingOption,
        extrapolation_option::ExtrapolationOption,
        buf::FitBuffer2D,
        S::Matrix{T},
        x1_start::T,
        x1_fin::T,
        x2_start::T,
        x2_fin::T;
        ϵ::T=eps(T) * 2,
    ) where {T<:AbstractFloat}

        A1, xs1 = create_query_cache(x1_start, x1_fin, size(S, 1), buf.buf_x1.N_padding)
        A2, xs2 = create_query_cache(x2_start, x2_fin, size(S, 2), buf.buf_x2.N_padding)

        c = zeros(T, get_coeffs_size(buf)) # allocates.
        _get_coeffs!(padding_option, extrapolation_option, c, buf, S, xs1, xs2, ϵ)

        return new{T}(c, A1, A2, x1_start, x1_fin, x2_start, x2_fin)
    end

    # has padding. Default to LinearPadding()
    function Interpolator2D(buf::FitBuffer2D, S::Matrix{T}, x1_start::T, x1_fin::T, x2_start::T, x2_fin::T; ϵ::T=eps(T) * 2) where {T<:AbstractFloat}
        return Interpolator2D(LinearPadding(), ConstantExtrapolation(), buf, S, x1_start, x1_fin, x2_start, x2_fin; ϵ=ϵ)
    end
end

# option is for dispatch. itp mutates, is output. Mutates buf, 
function update_itp!(
    padding_option::PaddingOption,
    extrapolation_option::ExtrapolationOption,
    itp::Interpolator2D,
    buf::FitBuffer2D,
    S::Matrix{T};
    ϵ::T=eps(T) * 2,
) where {T<:AbstractFloat}

    size(itp.coeffs) == get_coeffs_size(buf) || error("Size mismatch.")
    size(S, 2) == size(buf.S1, 2) || error("Size mismatch.")

    x1_start, x1_fin, x2_start, x2_fin = get_itp_interval(itp)
    xs1 = LinRange(x1_start, x1_fin, size(S, 1))
    xs2 = LinRange(x2_start, x2_fin, size(S, 2))
    _get_coeffs!(padding_option, extrapolation_option, itp.coeffs, buf, S, xs1, xs2, ϵ)

    return nothing
end

# convenince
function update_itp!(itp::Interpolator2D, buf::FitBuffer2D, S::Matrix{T}; ϵ::T=eps(T) * 2) where {T<:AbstractFloat}
    return update_itp!(LinearPadding(), ConstantExtrapolation(), itp, buf, S; ϵ=ϵ)
end


# X, buf1, buf2 mutates, are buffers.
# Y mutates, is output.
function _get_coeffs!(
    pading_option::PaddingOption,
    extrapolation_option::ExtrapolationOption,
    Y::AbstractMatrix{T},
    buf::FitBuffer2D,
    S::AbstractMatrix{T},
    xs1::LinRange,
    xs2::LinRange,
    ϵ::T,
) where {T<:AbstractFloat}

    size(Y, 1) == size(buf.S1, 1) || error("Size mismatch.")
    size(S, 2) == size(buf.S1, 2) || error("Size mismatch.")

    X = buf.S1
    for (xc, sc) in Iterators.zip(eachcol(X), eachcol(S))
        _get_coeffs!(pading_option, extrapolation_option, xc, buf.buf_x1, sc, xs1, ϵ)
    end

    for (yr, xr) in Iterators.zip(eachrow(Y), eachrow(X))

        _get_coeffs!(pading_option, extrapolation_option, yr, buf.buf_x2, xr, xs2, ϵ)
    end

    return nothing
end


function query2D(x1_in::T, x2_in::T, itp::Interpolator2D{T}) where {T<:AbstractFloat}
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

    # The 4x4 = 16 terms.
    out00 = C[begin+k1_lb, begin+k2_lb] * spline1_0 * spline2_0
    out10 = C[begin+k1_lb+1, begin+k2_lb] * spline1_1 * spline2_0
    out20 = C[begin+k1_lb+2, begin+k2_lb] * spline1_2 * spline2_0
    out30 = C[begin+k1_lb+3, begin+k2_lb] * spline1_3 * spline2_0

    out01 = C[begin+k1_lb, begin+k2_lb+1] * spline1_0 * spline2_1
    out11 = C[begin+k1_lb+1, begin+k2_lb+1] * spline1_1 * spline2_1
    out21 = C[begin+k1_lb+2, begin+k2_lb+1] * spline1_2 * spline2_1
    out31 = C[begin+k1_lb+3, begin+k2_lb+1] * spline1_3 * spline2_1

    out02 = C[begin+k1_lb, begin+k2_lb+2] * spline1_0 * spline2_2
    out12 = C[begin+k1_lb+1, begin+k2_lb+2] * spline1_1 * spline2_2
    out22 = C[begin+k1_lb+2, begin+k2_lb+2] * spline1_2 * spline2_2
    out32 = C[begin+k1_lb+3, begin+k2_lb+2] * spline1_3 * spline2_2

    out03 = C[begin+k1_lb, begin+k2_lb+3] * spline1_0 * spline2_3
    out13 = C[begin+k1_lb+1, begin+k2_lb+3] * spline1_1 * spline2_3
    out23 = C[begin+k1_lb+2, begin+k2_lb+3] * spline1_2 * spline2_3
    out33 = C[begin+k1_lb+3, begin+k2_lb+3] * spline1_3 * spline2_3

    return out00 + out10 + out20 + out30 +
           out01 + out11 + out21 + out31 +
           out02 + out12 + out22 + out32 +
           out03 + out13 + out23 + out33
end

