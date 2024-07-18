
# Box 2, Unser 1999.
function get_coeffs(s::Memory{T}, ϵ::T) where T <: AbstractFloat
    c_minus = Memory{T}(undef, length(s))
    _get_coeffs!(c_minus, s, ϵ)
    return c_minus
end

# Not safe for public use due to assumptions on `c_minus` and `s`.
# `c_mins` and `s` must be an AbstractVector with stride 1, 1-indexing.
function _get_coeffs!(c_minus, s, ϵ::T) where T <: AbstractFloat

    @assert length(c_minus) == length(s)
    @assert !isempty(s)

    z1 = T(sqrt(3)-2)
    
    k0 = clamp(round(Int, log(ϵ)/log(abs(z1))), 2, length(s))

    N = length(s)
    c_plus = Memory{T}(undef, N)
    fill!(c_plus, 0)
    fill!(c_minus, 0)
    
    # # Causal filtering to get intermediate coefficients.

    # Use the approximation by Unser to get this initial condition.
    c_plus[begin] = evalpoly(z1, view(s, 1:k0)) #/(1-z1^(2*N-2))
    
    for k = 2:N
        c_plus[k] = s[k] + z1*c_plus[k-1]
    end

    # # Acausal filtering to get final coefficients.
    
    # terminal condition.
    c_minus[end] = z1/(1-z1^2) * (c_plus[end] + z1*c_plus[end-1])

    for k = N-1:-1:1
        c_minus[k] = z1*(c_minus[k+1] - c_plus[k])
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
end

# Specialized for 1 <= abs(x) < 2
function eval_cubic_spline_in12(x::AbstractFloat)
    abs_x = abs(x)
    return ((2-abs_x)^3)/6
end

# For x ∈ [a,b], where (a,b,d) were used to create A.
# Full convolution.
function query(x_in::T, c::Memory{T}, A::IntervalConversion{T}) where T <: AbstractFloat
    x = to_std_interval(x_in, A.a, A.d_div_bma)

    out = zero(T)
    for k in eachindex(c)
        out += c[k]*eval_cubic_spline(x-k+1)
    end
    return out
end

# For x ∈ [a,b], where (a,b,d) were used to create A.
# Use only specialized splines.
function query_interior(x_in::T, c::Memory{T}, A::IntervalConversion{T}) where T <: AbstractFloat
    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    x = clamp(x, 2, length(c)-1-2)

    k_lb = ceil(Int, x-2)
    
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

# package up
struct Interpolator1D{T <: AbstractFloat}
    coeffs::Memory{T}
    query_cache::IntervalConversion{T}
end

function Interpolator1D(s::Memory{T}, x_start::T, x_fin::T; ϵ::T = eps(T)*100) where T <: AbstractFloat
    
    A = IntervalConversion(x_start, x_fin, length(s))
    c = get_coeffs(s, ϵ)

    return Interpolator1D(c, A)
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
    valid_lb = t[4]
    valid_ub = t[end-3]

    return valid_lb, valid_ub
end

function query_interior(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat
    return query_interior(x_in, itp.coeffs, itp.query_cache)
end

######### 2D

function get_coeffs(S::Matrix{T}, ϵ::T) where T <: AbstractFloat
    
    X = zeros(T, size(S))
    for (xc,sc) in Iterators.zip(eachcol(X),eachcol(S))
        _get_coeffs!(xc, sc, ϵ)
    end

    Y = zeros(T, size(X))
    for (yr,xr) in Iterators.zip(eachrow(Y),eachrow(X))
        _get_coeffs!(yr, xr, ϵ)
    end

    return reshape(Memory{T}(vec(Y)), size(Y))
end

# Note: could use a macro to generate the 25 terms.
# x1 is the row direction, x2 is the column direction.
function query_interior(x1_in::T, x2_in::T, C, A1::IntervalConversion{T}, A2::IntervalConversion{T}) where T <: AbstractFloat
    
    # # Transform clamp the input coordinates.
    x1 = to_std_interval(x1_in, A1.a, A1.d_div_bma)
    x2 = to_std_interval(x2_in, A2.a, A2.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    x1 = clamp(x1, 2, size(C,1)-1-2)
    x2 = clamp(x2, 2, size(C,2)-1-2)

    # # Compute the 4x4 = 16 terms.
    # explicitly write it all out here to encourage compiler & LLVM to optimize, as the instructions can be parallelized.
    # perhaps in the future, the compiler and LLVM can do optimization for sub-level of functions reliably.
    k1_lb = ceil(Int, x1-2)
    k2_lb = ceil(Int, x2-2)
    
    # ## Fix column
    k2 = k2_lb
    x2_minus_k2 = x2 - k2

    # ### iterate rows.
    k1 = k1_lb
    out11 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)

    k1 = k1_lb + 1
    out21 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)

    k1 = k1_lb + 2
    out31 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)
    
    k1 = k1_lb + 3
    out41 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)

    # ## Fix column
    k2 = k2_lb + 1
    x2_minus_k2 = x2 - k2

    # ### iterate rows.
    k1 = k1_lb
    out12 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)

    k1 = k1_lb + 1
    out22 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)

    k1 = k1_lb + 2
    out32 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)
    
    k1 = k1_lb + 3
    out42 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)

    # ## Fix column
    k2 = k2_lb + 2
    x2_minus_k2 = x2 - k2

    # ### iterate rows.
    k1 = k1_lb
    out13 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)

    k1 = k1_lb + 1
    out23 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)

    k1 = k1_lb + 2
    out33 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)
    
    k1 = k1_lb + 3
    out43 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in01(x2_minus_k2)
    
    # ## Fix column
    k2 = k2_lb + 3
    x2_minus_k2 = x2 - k2

    # ### iterate rows.
    k1 = k1_lb
    out14 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)

    k1 = k1_lb + 1
    out24 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)

    k1 = k1_lb + 2
    out34 = C[begin + k1,begin + k2]*eval_cubic_spline_in01(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)
    
    k1 = k1_lb + 3
    out44 = C[begin + k1,begin + k2]*eval_cubic_spline_in12(x1 - k1)*eval_cubic_spline_in12(x2_minus_k2)
    
    return out11 + out21 + out31 + out41 +
        out12 + out22 + out32 + out42 +
        out13 + out23 + out33 + out43 +
        out14 + out24 + out34 + out44
end

# package up.


struct Interpolator2D{T <: AbstractFloat, RT <: AbstractArray}
    coeffs::RT
    x1_query_cache::IntervalConversion{T}
    x2_query_cache::IntervalConversion{T}
end
# `RT` in v1.11, this is Base.ReshapedArray{Float64, 2, Memory{Float64}, Tuple{}}. Leave it as a parametric type in case Base changes this data type in the future.

function Interpolator2D(S::Matrix{T}, x1_start::T, x1_fin::T, x2_start::T, x2_fin::T; ϵ::T = eps(T)*100) where T <: AbstractFloat
    
    A1 = IntervalConversion(x1_start, x1_fin, size(S,1))
    A2 = IntervalConversion(x2_start, x2_fin, size(S,2))
    c = get_coeffs(S, ϵ)

    return Interpolator2D(c, A1, A2)
end

function get_query_interval(itp::Interpolator2D)
    x1_lb, x1_ub = get_query_interval(itp.x1_query_cache, size(itp.coeffs,1))
    x2_lb, x2_ub = get_query_interval(itp.x2_query_cache, size(itp.coeffs,2))
    return x1_lb, x1_ub, x2_lb, x2_ub
end

function query_interior(u1::T, u2::T, itp::Interpolator2D{T}) where T <: AbstractFloat
    return query_interior(u1, u2, itp.coeffs, itp.x1_query_cache, itp.x2_query_cache)
end

