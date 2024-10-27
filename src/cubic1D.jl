# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# abstract type AbstractInterpolator{T} end
# abstract type AbstractFitBuffer{T} end

# # convinence
# function update_itp!(itp::AbstractInterpolator, buf::AbstractFitBuffer, args...; ϵ = eps(T)*10) where T
#     return update_itp!(LinearPadding(), itp, buf, args...; ϵ = ϵ)
# end

abstract type AbstractInterpolator1D end
# interface requirements: have coeffs, x_start, x_fin as field names.
# have update_itp! as method

function get_itp_interval(itp::AbstractInterpolator1D)
    return itp.x_start, itp.x_fin
end

abstract type AbstractInterpolator2D end
# interface requirements: have coeffs, x1_start, x2_fin, x2_start, x2_fin, as field names.
# have update_itp! as method

function get_itp_interval(itp::AbstractInterpolator2D)
    return itp.x1_start, itp.x1_fin, itp.x2_start, itp.x2_fin
end

function update_coeffs!(itp::Union{AbstractInterpolator1D, AbstractInterpolator2D}, v::AbstractArray)
    copyto!(itp.coeffs, v)
    return nothing
end

function get_coeffs_len(itp::Union{AbstractInterpolator1D, AbstractInterpolator2D})
    return length(itp.coeffs)
end

# pre-allocated version
struct FitBuffer1D{T <: AbstractFloat}
    #c_minus::Memory{T}
    c_plus::Memory{T}
    s_extension::Memory{T}
    N_padding::Int

    function FitBuffer1D(::Type{T}, N::Integer; N_padding::Integer = 10) where T <: AbstractFloat

        N_padding >= 5 || error("N_padding must be larger or equal than 5. Two for the query system used in this library, three for transition to constant extrapolation.")
        N > 2 || error("N is the number of samples, and it must be larger than 2.")

        M = N + 2*N_padding # 2 extra samples at the each boundary.
        return new{T}(Memory{T}(undef, M), Memory{T}(undef, M), Int(N_padding))
    end
end

function get_data_length(s::FitBuffer1D)
    return length(s.s_extension) - s.N_padding*2 # padded two samples at each boundary.
end

function get_coeffs_length(s::FitBuffer1D)
    return length(s.c_plus)
end

# Box 2, Unser 1999.
function get_coeffs!(
    pad_option::PaddingOption,
    ex_option::ExtrapolationOption,
    c_minus::AbstractVector{T},
    buf::FitBuffer1D{T},
    s::AbstractVector{T},
    xs::LinRange,
    ϵ::T,
    ) where T <: AbstractFloat
    
    # @show get_data_length(buf), length(s)
    get_data_length(buf) == length(s) || error("Length mismatch.")
    length(s) > 2 || error("Must provide at least 2 input samples.")

    s_ext, c_plus = buf.s_extension, buf.c_plus

    # debug
    fill!(s_ext, 0)

    # pad s_extension.
    Np = buf.N_padding
    st_ind = Np + 1
    fin_ind = length(s_ext) - Np
    v = view(s_ext, st_ind:fin_ind)
    copy!(v, s)

    _pad_samples!(pad_option, s_ext, s, Np, xs)

    # get coefficient on the extended sample set.
    #_get_coeffs!(c_minus, buf.c_plus, s_ext, ϵ)
    _get_coeffs!(c_minus, c_plus, s_ext, ϵ)

    _post_process_coeffs!(ex_option, c_minus, Np)
    return nothing
end

# c must be 1-indexing, stride 1
function _post_process_coeffs!(::ConstantExtrapolation, c::AbstractVector, Np::Integer)

    M0 = Np-1 # ensure the query matches the samples.
    v = view(c, 1:M0)
    tmp = c[M0]
    #tmp = sum(v)/length(v)
    fill!(v, tmp)
    
    v = view(c, (length(c)-M0+1):length(c))
    tmp = c[(length(c)-M0+1)]
    #tmp = sum(v)/length(v)
    fill!(v, tmp)

    return nothing
end

# c must be 1-indexing, stride 1
function _post_process_coeffs!(::ZeroExtrapolation, c::AbstractVector, Np::Integer)

    M0 = Np-1
    v = view(c, 1:M0)
    fill!(v, 0)
    
    v = view(c, (length(c)-M0+1):length(c))
    fill!(v, 0)

    return nothing
end

# Box 2, Unser 1999.
function get_coeffs!(c_minus::Union{Memory{T},Vector{T}}, s::Union{Memory{T},Vector{T}}, ϵ::T) where T <: AbstractFloat
    length(c_minus) == length(s) || error("Length mismatch.")
    
    c_plus = Memory{T}(undef, length(s))
    _get_coeffs!(c_minus, c_plus, s, ϵ)
    return nothing
end

# Not safe for public use due to assumptions on `c_minus` and `s`.
# `c_mins` and `s` must be an AbstractVector with stride 1, 1-indexing.
# c_minus and s should be 1-indexing with stride 1.
function _get_coeffs!(
    c_minus::AbstractVector{T}, # output, mutates
    c_plus::Memory{T}, # buffer, mutates.
    s::AbstractVector{T}, 
    ϵ::T,
    ) where T <: AbstractFloat

    length(c_minus) == length(s) || error("Length mismatch.")
    !isempty(s) || error("Empty samples array.")

    z1 = sqrt(T(3)) - 2
    
    k0 = clamp(round(Int, log(ϵ)/log(abs(z1))), 2, length(s)-1)

    N = length(s)
    #c_plus = Memory{T}(undef, N) # TODO allocates.
    fill!(c_plus, 0)
    fill!(c_minus, 0)
    
    # # Causal filtering to get intermediate coefficients.

    # Use the approximation by Unser to get this initial condition.
    c_plus[begin] = evalpoly(z1, view(s, 1:k0+1)) /(1-z1^(2*N-2))
    
    # for k = 2:N
    #     c_plus[k] = s[k] + z1*c_plus[k-1]
    # end
    for k = 1:N-1
        c_plus[begin+k] = s[begin+k] + z1*c_plus[begin+k-1]
    end

    # # Acausal filtering to get final coefficients.
    
    # terminal condition.
    c_minus[end] = z1/(1-z1^2) * (c_plus[end] + z1*c_plus[end-1])

    # for k = N-1:-1:1
    #     c_minus[k] = z1*(c_minus[k+1] - c_plus[k])
    # end
    for k = N-2:-1:0
        c_minus[begin+k] = z1*(c_minus[begin+k+1] - c_plus[begin+k])
    end

    # go from c_minus to c.
    for i in eachindex(c_minus)
        c_minus[i] *= 6
    end

    return nothing
end

# # Eq 1, Unser 1999. Only for x ∈ [0, N-1].
# function query(x::T, c::Memory{T}) where T <: AbstractFloat

#     out = zero(T)
#     for k in 0:(length(c)-1)
#         out += c[begin+k]*eval_cubic_spline(x-k)
#     end

#     return out
# end

# Eq 6, Unser 1999. This is only used for reference and testing.
function eval_cubic_spline(x::T) where T <: AbstractFloat
    
    abs_x = abs(x)
    
    if 0 <= abs_x <1
        return 2/3 - abs_x^2 + (abs_x^3)/2
    
    elseif 1 <= abs_x < 2
        return ((2-abs_x)^3)/6
    end

    return zero(T)
end

# Specialized for 0 <= abs(x) < 1
function eval_cubic_spline_in01(x::T) where T <: AbstractFloat
    abs_x = abs(x)
    return twothirds(T) - abs_x^2 + (abs_x^3)/2

    #abs_x_sq = abs_x*abs_x
    #return twothirds(T) - abs_x_sq + (abs_x_sq*abs_x)/2
end

# Specialized for 1 <= abs(x) < 2
function eval_cubic_spline_in12(x::AbstractFloat)
    abs_x = abs(x)
    return ((2-abs_x)^3)/6
end

function create_query_cache(x_start, x_fin, N, Np)

    tmp_r = LinRange(x_start, x_fin, N)
    Δr = step(tmp_r)
    N_ext = N + 2*Np
    x_range_ext = LinRange(x_start - Δr*Np, x_fin + Δr*Np, N_ext)

    A = IntervalConversion(first(x_range_ext), last(x_range_ext), N_ext)
    return A, tmp_r
end

# package up
struct Interpolator1D{T <: AbstractFloat} <: AbstractInterpolator1D
    coeffs::Memory{T}
    query_cache::IntervalConversion{T}
    
    x_start::T
    x_fin::T

    # # The boundary 2 samples (e.g. in the Δx*2 region).won't match the corresponding boundary 2 samples in s.
    # function Interpolator1D(s::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
    #     A = IntervalConversion(x_start, x_fin, length(s))
    #     c = Memory{T}(undef, length(s))
    #     get_coeffs!(c, s, ϵ)
    #     return new{T}(c, A, x_start, x_fin)
    # end

    # Mutates buf, option is for dispatch.
    function Interpolator1D(
        pading_option::PaddingOption,
        extrapolation_option::ExtrapolationOption,
        buf::FitBuffer1D, s::Union{Memory{T},Vector{T}},
        x_start::T,
        x_fin::T;
        ϵ::T = eps(T)*2,
        ) where T <: AbstractFloat

        x_start < x_fin || error("x_start must be strictly smaller than x_fin.")
        
        Np = buf.N_padding

        # N = length(s)
        # tmp_r = LinRange(x_start, x_fin, N)
        # Δr = step(tmp_r)
        # N_ext = N + 2*Np
        # x_range_ext = LinRange(x_start - Δr*Np, x_fin + Δr*Np, N_ext)
        # A = IntervalConversion(first(x_range_ext), last(x_range_ext), N_ext)
        A, tmp_r = create_query_cache(x_start, x_fin, length(s), Np)
        
        c = Memory{T}(undef, get_coeffs_length(buf))
        get_coeffs!(pading_option, extrapolation_option, c, buf, s, tmp_r, ϵ)

        return new{T}(c, A, x_start, x_fin)
    end

    # convinence constructor. Mutates buf
    function Interpolator1D(buf::FitBuffer1D, s::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*2) where T <: AbstractFloat
        return Interpolator1D(LinearPadding(), ConstantExtrapolation(), buf, s, x_start, x_fin; ϵ = ϵ)
    end    
end

# option is for dispatch. itp mutates, is output. Mutates buf, 
function update_itp!(padding_option::PaddingOption, extrapolation_option::ExtrapolationOption, itp::Interpolator1D, buf::FitBuffer1D, s::Union{Memory{T},Vector{T}}; ϵ::T = eps(T)*2) where T <: AbstractFloat
    length(s) == get_data_length(buf) || error("Length mismatch.")

    x_start, x_fin = get_itp_interval(itp)
    get_coeffs!(padding_option, extrapolation_option, itp.coeffs, buf, s, LinRange(x_start, x_fin, length(s)), ϵ)
    return nothing
end

# convenince
function update_itp!(itp::Interpolator1D, buf::FitBuffer1D, s::Union{Memory{T},Vector{T}}; ϵ::T = eps(T)*2) where T <: AbstractFloat
    return update_itp!(LinearPadding(), ConstantExtrapolation(), itp, buf, s; ϵ = ϵ)
end


function query1D(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat
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
    out1 = c[begin + k]*eval_cubic_spline_in12(x - k)

    k = k_lb + 1
    out2 = c[begin + k]*eval_cubic_spline_in01(x - k)

    k = k_lb + 2
    out3 = c[begin + k]*eval_cubic_spline_in01(x - k)
    
    k = k_lb + 3
    out4 = c[begin + k]*eval_cubic_spline_in12(x - k)

    return out1 + out2 + out3 + out4
end
