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
