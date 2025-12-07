# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# minimalist code for an interpolation-only spline surrogate. Does not fit the coefficients, which are all initialized to zero.


PLT.close("all")
fig_num = 1

T = Float64
rng = Random.Xoshiro(0)

num_coeffs = 100

#N_obs = 25
N_obs = num_coeffs
Nq = 10_000

lb = T(-10)
ub = T(10)

Δx = (ub - lb) / (N_obs - 1)

X = LinRange(lb, ub, N_obs)
f = xx -> (sinc(xx) + 1)
f_X = f.(X)

itp = ITP.Interpolator1D(lb, ub, num_coeffs)
itp.coeffs .= f.(ITP.get_coeffs_range(itp)) # rough estimate for coefficients.

q = tt -> ITP.itp1D(tt, itp)
dq = tt -> ITP.itp1D_derivative1(tt, itp)
dq_FD = tt -> FiniteDiff.finite_difference_derivative(q, tt)

query_lb, query_ub = ITP.get_query_bounds(itp)
@assert abs(query_lb - lb) < eps(T)
@assert abs(query_ub - ub) < eps(T)

xqs = LinRange(query_lb, query_ub, Nq)
q_xqs = q.(xqs)
dq_xqs = dq.(xqs)
dq_FD_xqs = dq_FD.(xqs)

@show norm(dq_FD_xqs - dq_xqs) / norm(dq_FD_xqs)

PLT.figure(fig_num)
fig_num += 1
PLT.plot(X, f_X, label = "data", "x")
PLT.plot(xqs, q_xqs, label = "itp")
PLT.legend()
PLT.title("spline surrogate")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(X, dq.(X), label = "data positions", "s")
PLT.plot(xqs, dq_xqs, label = "itp derivative", linewidth = 3)
PLT.plot(xqs, dq_FD_xqs, "--", label = "FiniteDiff.jl", linewidth = 3)
PLT.legend()
PLT.title("derivative")

println("Make sure the boundaries are valid")
@show ITP.get_query_bounds(itp)
@show q(lb), q(ub), dq(lb), dq(ub)

@assert isfinite(q(lb))
@assert isfinite(dq(ub))
@assert isfinite(q(lb))
@assert isfinite(dq(ub))

# # Parameter gradient
# Treat the coefficients in `itp` as a set of parameters, and get the gradient wrt those parameters for the function `itp1D`.

derivative_rtol = T(1.0e-9)
num_test_positions = 4321
num_test_per_position = 100

dq_c = Memory{T}(zeros(T, ITP.get_num_coeffs(itp)))

function eval_fdf!(grad::AbstractVector, itp::ITP.Interpolator1D, c::AbstractVector, x0::AbstractFloat)
    length(grad) == length(c) || error("Length mismatch.")

    ITP.update_coeffs!(itp, c)
    q_x0, ind1, ind2, ind3, ind4, grad1, grad2, grad3, grad4 = ITP.itp1D_parameter_derivatives(x0, itp)

    fill!(grad, 0)
    grad[ind1] = grad1
    grad[ind2] = grad2
    grad[ind3] = grad3
    grad[ind4] = grad4

    return q_x0
end

function eval_f!(itp::ITP.Interpolator1D, c::AbstractVector, x0::AbstractFloat)
    ITP.update_coeffs!(itp, c)
    return ITP.itp1D(x0, itp)
end

x0_set = LinRange(lb, ub, num_test_positions)

for x0 in x0_set

    for _ in 1:num_test_per_position
        c0 = randn(rng, T, num_coeffs)

        qc = xx -> eval_f!(itp, xx, x0)
        dq_c_FD = FiniteDiff.finite_difference_gradient(qc, c0)

        q_c0 = eval_fdf!(dq_c, itp, c0, x0)

        @assert norm(dq_c_FD - dq_c) / norm(dq_c_FD) < derivative_rtol
    end
end

nothing
