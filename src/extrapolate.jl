# non-constant extrapolation, via function joining.

struct QuadraticCallable{T<:AbstractFloat}
    a::T
    b::T
    c::T

    function QuadraticCallable(itp::Interpolator1D, x0::T) where {T<:AbstractFloat}
        F0 = query1D(x0, itp)
        F1 = query1D_derivative1(x0, itp)
        F2 = query1D_derivative2(x0, itp)
        a, b, c = get_parabola(x0, F0, F1, F2)
        return new{T}(a, b, c)
    end
end
function (Q::QuadraticCallable)(x::Real)
    return Q.a * x^2 + Q.b * x + Q.c
end

function get_parabola(x0, F0, F1, F2)
    a = F2 / 2
    b = F1 - 2 * a * x0
    c = F0 - a * x0^2 - b * x0
    return a, b, c
end

struct ExtrapolatorCore1D{T<:AbstractFloat}

    # common to all extrapolators.
    itp::Interpolator1D{T}

    lb_st::T
    lb_fin::T

    ub_st::T
    ub_fin::T

    function ExtrapolatorCore1D(itp::Interpolator1D{T}; num_transition_samples::Integer=1) where {T<:AbstractFloat}
        num_transition_samples > 0 || error("num_transition_samples must be a positive integer.")

        A = itp.query_cache

        # construct the transition regions. One region for the lower boundary (lb), one for the upper boundary (ub).
        lb_fin, ub_st = get_itp_interval(itp)

        # float_ind_* means floating-point, but represents an index/integer.
        float_ind_lb_fin = to_std_interval(lb_fin, A.a, A.d_div_bma)
        float_ind_ub_st = to_std_interval(ub_st, A.a, A.d_div_bma)

        float_ind_lb_st = float_ind_lb_fin - num_transition_samples
        float_ind_ub_fin = float_ind_ub_st + num_transition_samples
        lb_st = from_std_interval(float_ind_lb_st, A.a, A.d_div_bma)
        ub_fin = from_std_interval(float_ind_ub_fin, A.a, A.d_div_bma)

        return new{T}(itp, lb_st, lb_fin, ub_st, ub_fin)
    end
end

struct QuadraticExtrapolator1D{T<:AbstractFloat}
    core::ExtrapolatorCore1D{T}

    # the lb and ub extrapolation functions.
    lb::QuadraticCallable{T}
    ub::QuadraticCallable{T}

    function QuadraticExtrapolator1D(itp::Interpolator1D{T}; num_transition_samples::Integer=1) where {T<:AbstractFloat}
        num_transition_samples > 0 || error("num_transition_samples must be positive.")
        core = ExtrapolatorCore1D(itp; num_transition_samples)

        # match 0-th, 1-st, 2-nd derivatives at interior query boundaries.
        lb = QuadraticCallable(itp, core.lb_fin)
        ub = QuadraticCallable(itp, core.ub_st)
        return new{T}(core, lb, ub)
    end
end

# # TODO placeholder for now.
# function eval_lb_ex_func(x)
#     x
# end

# # TODO placeholder for now.
# # fit parabola with position and derivatives at transition start.
# # This makes sure query1D and extrapolation func intersects at transition start time.
# function eval_ub_ex_func(x, a, b, c)
#     #return 1.98
#     return a*x^2 + b*x + c
# end

function query1D(x_in::T, etp::QuadraticExtrapolator1D{T}) where {T<:AbstractFloat}

    itp = etp.core.itp
    b, c = get_itp_interval(itp)
    # assume this is run most often.
    if b < x_in < c
        return query1D(x_in, itp)
    end

    core, lb_callable, ub_callable = etp.core, etp.lb, etp.ub
    a, d = core.lb_st, core.ub_fin

    # # Join function algorithm: https://www.youtube.com/watch?v=vD5g8aVscUI&t=523s

    if d < x_in # upper boundary extrapolation.
        return ub_callable(x_in)

    elseif c < x_in # upper boundary transition.

        f_x = query1D(x_in, itp)
        g_x = ub_callable(x_in)

        return _eval_transition(x_in, c, d, f_x, g_x)

    elseif a < x_in # lower boundary transition.

        g_x = query1D(x_in, itp)
        f_x = lb_callable(x_in)
        return _eval_transition(x_in, a, b, f_x, g_x)
    end

    return lb_callable(x_in) # lower boundary extrapolation.
end

function _eval_transition(x, c, d, f_x, g_x)
    # assumes c <= x <= d
    w_x = compute_transition_weight(x, c, d)

    return (1 - w_x) * f_x + w_x * g_x
end


function compute_transition_weight(x::T, a::T, b::T) where {T<:AbstractFloat}
    return compute_transition_weight((x - a) / (x - b))
end

function compute_transition_weight(x::T) where {T<:AbstractFloat}
    if x < 0
        return zero(T)
    end

    if x < 1
        p_x = eval_psi(x)
        return p_x / (p_x + eval_psi(one(T) - x))
    end

    return one(T)
end

function eval_psi(x::T) where {T<:AbstractFloat}
    if x > 0
        return exp(-1 / x)
    end
    return zero(T)
end
