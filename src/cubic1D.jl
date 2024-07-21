abstract type AbstractInterpolator1D{T} end
abstract type AbstractInterpolator2D{T} end

# Box 2, Unser 1999.
function get_coeffs(s::Union{Memory{T},Vector{T}}, ϵ::T) where T <: AbstractFloat
    c_minus = Memory{T}(undef, length(s))
    _get_coeffs!(c_minus, s, ϵ)
    return c_minus
end

# Not safe for public use due to assumptions on `c_minus` and `s`.
# `c_mins` and `s` must be an AbstractVector with stride 1, 1-indexing.
function _get_coeffs!(c_minus, s, ϵ::T) where T <: AbstractFloat

    length(c_minus) == length(s) || error("Length mismatch.")
    !isempty(s) || error("Empty samples array.")

    z1 = T(sqrt(3)-2)
    
    k0 = clamp(round(Int, log(ϵ)/log(abs(z1))), 2, length(s)-1)

    N = length(s)
    c_plus = Memory{T}(undef, N)
    fill!(c_plus, 0)
    fill!(c_minus, 0)
    
    # # Causal filtering to get intermediate coefficients.

    # Use the approximation by Unser to get this initial condition.
    c_plus[begin] = evalpoly(z1, view(s, 1:k0+1)) #/(1-z1^(2*N-2))
    
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

# Eq 1, Unser 1999. Only for x ∈ [0, N-1].
function query(x::T, c::Memory{T}) where T <: AbstractFloat

    out = zero(T)
    for k in 0:(length(c)-1)
        out += c[begin+k]*eval_cubic_spline(x-k)
    end

    return out
end

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

# For x ∈ [a,b], where (a,b,d) were used to create A.
# Full convolution. This is the reference formula, but we exclude the border 8 samples for query1D and query2D.
function query(x_in::T, c::Memory{T}, A::IntervalConversion{T}) where T <: AbstractFloat
    x = to_std_interval(x_in, A.a, A.d_div_bma)

    out = zero(T)
    for k in eachindex(c)
        out += c[k]*eval_cubic_spline(x-k+1)
    end
    return out
end


# package up
struct Interpolator1D{T <: AbstractFloat} <: AbstractInterpolator1D{T}
    coeffs::Memory{T}
    query_cache::IntervalConversion{T}

    function Interpolator1D(s::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*100) where T <: AbstractFloat
    
        A = IntervalConversion(x_start, x_fin, length(s))
        c = get_coeffs(s, ϵ)
    
        return new{T}(c, A)
    end
end

function get_query_interval(itp::Interpolator1D)
    return get_query_interval(itp.query_cache, length(itp.coeffs))
end

function get_query_interval(A::IntervalConversion, N)
    x_start = A.a # a
    
    d = N-1
    # d/(b-a) = d_div_bma
    x_fin = d/A.d_div_bma + x_start # b

    t = LinRange(x_start, x_fin, N)
    valid_lb = t[8]
    valid_ub = t[end-7]

    return valid_lb, valid_ub
end

# For x ∈ [a,b], where (a,b,d) were used to create A.
# Use only specialized splines.
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


# # Complex-valued 1D

struct Interpolator1DComplex{T <: AbstractFloat} <: AbstractInterpolator1D{T}
    
    real_coeffs::Memory{T}
    imag_coeffs::Memory{T}
    query_cache::IntervalConversion{T}

    function Interpolator1DComplex(s_real::Union{Memory{T},Vector{T}}, s_imag::Union{Memory{T},Vector{T}}, x_start::T, x_fin::T; ϵ::T = eps(T)*100) where T <: AbstractFloat
        size(s_real) == size(s_imag) || error("Size mismatch between the real and imaginary sample matrices.")
        
        A = IntervalConversion(x_start, x_fin, length(s_real))
        c_real = get_coeffs(s_real, ϵ)
        c_imag = get_coeffs(s_imag, ϵ)

        # sanity check. The following should pass since S_real and S_imag have the same size.
        # Keep it here to explicitly define invariants on real_coeffs and image_coeffs.
        size(c_real) == size(c_imag) || error("Size mismatch for coefficient matrix. Please file bug report.")
        return new{T}(c_real, c_imag, A)
    end
end

function get_query_interval(itp::Interpolator1DComplex)
    return get_query_interval(itp.query_cache, length(itp.real_coeffs))
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
