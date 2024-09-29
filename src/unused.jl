# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

function query_interior_old(x1_in::T, x2_in::T, C, A1::IntervalConversion{T}, A2::IntervalConversion{T}) where T <: AbstractFloat
    
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


# somehow, slower. I suspect it is the bounds check on C and the indices `k1_lb_1`, `k1_lb_2`, being determined at run-time, causing many different bounds check on `C`.
function query_interior2(x1_in::T, x2_in::T, C, A1::IntervalConversion{T}, A2::IntervalConversion{T}) where T <: AbstractFloat
    
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

    k1_lb_1 = k1_lb + 1
    k1_lb_2 = k1_lb + 2
    k1_lb_3 = k1_lb + 3
    
    k2_lb_1 = k2_lb + 1
    k2_lb_2 = k2_lb + 2
    k2_lb_3 = k2_lb + 3

    # Compute shared terms spline1_2 means direction x1, shift by 2 from the lower bound of the sampling lattice.
    # the point (x1, x2) is near the center of the 4x4 window on the sampling lattice, so shift 2 and 3 use the 01 version of the spline, similarly for shifts 0 and 4 for the 12 version.
    spline2_0 = eval_cubic_spline_in12(x2 - k2_lb)
    spline2_1 = eval_cubic_spline_in01(x2 - k2_lb_1)
    spline2_2 = eval_cubic_spline_in01(x2 - k2_lb_2)
    spline2_3 = eval_cubic_spline_in12(x2 - k2_lb_3)

    spline1_0 = eval_cubic_spline_in12(x1 - k1_lb)
    spline1_1 = eval_cubic_spline_in01(x1 - k1_lb_1)
    spline1_2 = eval_cubic_spline_in01(x1 - k1_lb_2)
    spline1_3 = eval_cubic_spline_in12(x1 - k1_lb_3)

    # The 4x4 = 16 terms.
    out00 = C[begin + k1_lb, begin + k2_lb]*spline1_0*spline2_0
    out10 = C[begin + k1_lb_1, begin + k2_lb]*spline1_1*spline2_0
    out20 = C[begin + k1_lb_2, begin + k2_lb]*spline1_2*spline2_0
    out30 = C[begin + k1_lb_3, begin + k2_lb]*spline1_3*spline2_0

    out01 = C[begin + k1_lb, begin + k2_lb_1]*spline1_0*spline2_1
    out11 = C[begin + k1_lb_1, begin + k2_lb_1]*spline1_1*spline2_1
    out21 = C[begin + k1_lb_2, begin + k2_lb_1]*spline1_2*spline2_1
    out31 = C[begin + k1_lb_3, begin + k2_lb_1]*spline1_3*spline2_1

    out02= C[begin + k1_lb, begin + k2_lb_2]*spline1_0*spline2_2
    out12 = C[begin + k1_lb_1, begin + k2_lb_2]*spline1_1*spline2_2
    out22 = C[begin + k1_lb_2, begin + k2_lb_2]*spline1_2*spline2_2
    out32 = C[begin + k1_lb_3, begin + k2_lb_2]*spline1_3*spline2_2

    out03= C[begin + k1_lb, begin + k2_lb_3]*spline1_0*spline2_3
    out13 = C[begin + k1_lb_1, begin + k2_lb_3]*spline1_1*spline2_3
    out23 = C[begin + k1_lb_2, begin + k2_lb_3]*spline1_2*spline2_3
    out33 = C[begin + k1_lb_3, begin + k2_lb_3]*spline1_3*spline2_3

    return out00 + out10 + out20 + out30 +
        out01 + out11 + out21 + out31 +
        out02 + out12 + out22 + out32 +
        out03 + out13 + out23 + out33
end


function query_interior2(u1::T, u2::T, itp::Interpolator2D{T}) where T <: AbstractFloat
    return query_interior2(u1, u2, itp.coeffs, itp.x1_query_cache, itp.x2_query_cache)
end

function query_interior_original(u1::T, u2::T, itp::Interpolator2D{T}) where T <: AbstractFloat
    return query_interior_original(u1, u2, itp.coeffs, itp.x1_query_cache, itp.x2_query_cache)
end

# Note: could use a macro to generate the 25 terms.
# x1 is the row direction, x2 is the column direction.
function query_interior_original(x1_in::T, x2_in::T, C, A1::IntervalConversion{T}, A2::IntervalConversion{T}) where T <: AbstractFloat
    
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
    k2_lb = ceil(Int, x2-2) # takes a long time.

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
    out00 = C[begin + k1_lb, begin + k2_lb]*spline1_0*spline2_0
    out10 = C[begin + k1_lb + 1, begin + k2_lb]*spline1_1*spline2_0
    out20 = C[begin + k1_lb + 2, begin + k2_lb]*spline1_2*spline2_0
    out30 = C[begin + k1_lb + 3, begin + k2_lb]*spline1_3*spline2_0

    out01 = C[begin + k1_lb, begin + k2_lb + 1]*spline1_0*spline2_1
    out11 = C[begin + k1_lb + 1, begin + k2_lb + 1]*spline1_1*spline2_1
    out21 = C[begin + k1_lb + 2, begin + k2_lb + 1]*spline1_2*spline2_1
    out31 = C[begin + k1_lb + 3, begin + k2_lb + 1]*spline1_3*spline2_1

    out02= C[begin + k1_lb, begin + k2_lb + 2]*spline1_0*spline2_2
    out12 = C[begin + k1_lb + 1, begin + k2_lb + 2]*spline1_1*spline2_2
    out22 = C[begin + k1_lb + 2, begin + k2_lb + 2]*spline1_2*spline2_2
    out32 = C[begin + k1_lb + 3, begin + k2_lb + 2]*spline1_3*spline2_2

    out03= C[begin + k1_lb, begin + k2_lb + 3]*spline1_0*spline2_3
    out13 = C[begin + k1_lb + 1, begin + k2_lb + 3]*spline1_1*spline2_3
    out23 = C[begin + k1_lb + 2, begin + k2_lb + 3]*spline1_2*spline2_3
    out33 = C[begin + k1_lb + 3, begin + k2_lb + 3]*spline1_3*spline2_3

    return out00 + out10 + out20 + out30 +
        out01 + out11 + out21 + out31 +
        out02 + out12 + out22 + out32 +
        out03 + out13 + out23 + out33
end


# # region extensions.

# assumes a < b
function transition_weight(x, a, b)
    return transition_weight((x-a)/(b-a))
end

function eval_ψ(x::T) where T <: AbstractFloat
    if x > 0
        return exp(-one(T)/x)
    end
    return zero(T)
end

function transition_weight(x::T) where T <: AbstractFloat

    if x < 0
        return zero(T)
    end

    if x > 1
        return one(T)
    end

    ψ_x = eval_ψ(x)
    return ψ_x/(ψ_x + eval_ψ(one(T)-x))
end

# no longer used. no mirror image required.
#### move this idea to another repo.
function eval_reflection_etp(ts::LinRange, samples::AbstractVector{T}) where T <: AbstractFloat

    Δt = step(ts)
    order = 4

    # # the larger boundary
    xs = ts[length(ts)-4:length(ts)]
    ys = view(samples, length(ts)-4:length(ts))
    w0, w1, w2, w3, w4 = setup_lagrange4(xs)
    
    a = xs[end]
    x1 = a + Δt
    x2 = a + 2*Δt
    x3 = a + 3*Δt
    b = a + 4*Δt

    f_x1 = lagrange4_etp(x1, xs, ys, w0, w1, w2, w3, w4)
    g_x1 = lagrange4_etp(x3, xs, ys, w0, w1, w2, w3, w4)
    ϕ_x1 = transition_weight(x1, a, b)
    ext_large_x1 = (one(T) - ϕ_x1)*f_x1 + ϕ_x1*g_x1
    
    f_x2 = lagrange4_etp(x2, xs, ys, w0, w1, w2, w3, w4)
    g_x2 = lagrange4_etp(x2, xs, ys, w0, w1, w2, w3, w4)
    ϕ_x2 = transition_weight(x2, a, b)
    ext_large_x2 = (one(T) - ϕ_x2)*f_x2 + ϕ_x2*g_x2

    # f(x) # x < a
    # (1-ϕ(x))*f(x) + ϕ(x)*g(x)
    # g(x) # x > b

    # # the larger boundary
    xs = ts[1:(order+1)]
    ys = view(samples, 1:(order+1))
    w0, w1, w2, w3, w4 = setup_lagrange4(xs)
    
    b = xs[begin]
    x3 = b - Δt
    x2 = b - 2*Δt
    x1 = b - 3*Δt
    a = b - 4*Δt

    f_x3 = lagrange4_etp(x3, xs, ys, w0, w1, w2, w3, w4)
    g_x3 = lagrange4_etp(x1, xs, ys, w0, w1, w2, w3, w4)
    ϕ_x3 = transition_weight(x3, a, b)
    ext_small_x3 = (one(T) - ϕ_x3)*f_x3 + ϕ_x3*g_x3
    
    f_x2 = lagrange4_etp(x2, xs, ys, w0, w1, w2, w3, w4)
    g_x2 = lagrange4_etp(x2, xs, ys, w0, w1, w2, w3, w4)
    ϕ_x2 = transition_weight(x2, a, b)
    ext_small_x2 = (one(T) - ϕ_x2)*f_x2 + ϕ_x2*g_x2

    return ext_small_x3, ext_small_x2, ext_large_x1, ext_large_x2
end


#### assuming coefficients are zero outside of signal boundary.

# For x ∈ [a,b], where (a,b,d) were used to create A.
# Use only specialized splines.
function query1D_formula(x_in::T, itp::Interpolator1D{T}) where T <: AbstractFloat
    c, A = itp.coeffs, itp.query_cache

    x = to_std_interval(x_in, A.a, A.d_div_bma)
    
    # clamp x to be within bounds.
    # This has the effect that the border (2 pixels, inclusive) and extrapolation behavior is clamped (i.e. a constant) to the closest non-border interpolation.
    # this ensures no out-of-bounds access to `c` for index `kp1`.
    #x = clamp(x, 2, length(c)-1-2)
    # x_lb = T(2)
    
    # x_ub = T(size(c,1)-1-2)
    # if x < x_lb
    #     x = x_lb
    # elseif x > x_ub
    #     x = x_ub
    # end

    #k_lb = ceil(Int, x-2)
    k_lb = trunc(Int, x)-1
    
    k = k_lb
    #out1 = c[begin + k]*eval_cubic_spline(x - k) # non-specialized.
    p = zero(T)
    if 1 <= 1+k <= length(c)
        p =c[begin + k]
    end 
    out1 = p*eval_cubic_spline_in12(x - k)

    k = k_lb + 1
    p = zero(T)
    if 1 <= 1+k <= length(c)
        p =c[begin + k]
    end 
    out2 = p*eval_cubic_spline_in01(x - k)

    k = k_lb + 2
    p = zero(T)
    if 1 <= 1+k <= length(c)
        p =c[begin + k]
    end 
    out3 = p*eval_cubic_spline_in01(x - k)
    
    k = k_lb + 3
    p = zero(T)
    if 1 <= 1+k <= length(c)
        p =c[begin + k]
    end 
    out4 = p*eval_cubic_spline_in12(x - k)

    return out1 + out2 + out3 + out4
end


#### backup

function get_coeffs!(Y::Matrix{T}, X::Matrix, buf1::AbstractVector, buf2::AbstractVector, S::Matrix{T}, ϵ::T) where T <: AbstractFloat
    length(buf1) == size(S,1) || error("Length mismatch.")
    length(buf2) == size(S,2) || error("Length mismatch.")
    size(Y) == size(X) == size(S) || error("Size mismatch.")

    X = zeros(T, size(S))
    for (xc,sc) in Iterators.zip(eachcol(X),eachcol(S))
        _get_coeffs!(xc, buf1, sc, ϵ)
    end

    Y = zeros(T, size(X))
    for (yr,xr) in Iterators.zip(eachrow(Y),eachrow(X))
        _get_coeffs!(yr, buf2, xr, ϵ)
    end

    #return reshape(Memory{T}(vec(Y)), size(Y))
    return nothing
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


#### derivatives
# These are incorrect.

# specialized for 1 < abs(x) < 2 # incorrect result. possibly due to our choice of query bounds.
function eval_d2B3_in12(x::T) where T <: AbstractFloat
    if x > 0
        return -(2-x)
        #return 2 + x
    end

    return 2 + x
    #return -(2-x)
end

# specialized for 1 < abs(x) < 2 # incorrect result. possibly due to our choice of query bounds.
# specialized for 0 <= abs(x) < 1
function eval_d2B3_in01(x::T) where T <: AbstractFloat
    if x > 0
        return -2 + 3*x
        #return -2 - 3*x
    end

    return -2 - 3*x
    #return -2 + 3*x
end

# suspected incorrect
# specialized for 1 < abs(x) < 2
function eval_dB3_in12(x::T) where T <: AbstractFloat
    if x > 0
        return -half(T)*(2-x)^2
    end

    return half(T)*(2+x)^2
end

# suspected incorrect
# specialized for 0 <= abs(x) < 1
function eval_dB3_in01(x::T) where T <: AbstractFloat
    if x > 0
        return -2*x + T(3)/T(2) * x^2
    end

    return -2*x -T(3)/T(2) *x^2
end