# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>


if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

rng = Random.Xoshiro(0)

const T = Float64
#const T = Float32

#ϵ = T(1e-16)
ϵ = T(1.0e-8)
Nq = 757
N = 1000
#N = 7 # debug
a = T(-10)
b = T(3.4)
t_range = LinRange(a, b, N)
#f = xx->sinc(sin((xx/5)^3+(xx/20)^2+4.32)) # not bad
#f = xx->sin(xx) # not bad.
#f = xx->(exp(-(1/17)*(xx-3)^2) +55) # runge effect near right boundary.
f = xx -> exp(-(1 / 17) * (xx - 3)^2) # runge effect near right boundary.
#f = xx->(exp(-(1/17)*(xx-3)^2) + sin(xx-1.2)) # runge effect near right boundary.

#f = xx->exp(-(1/17)*(xx/3)^2) # not bad.
#f = xx->sinc(xx)

A = ITP.IntervalConversion(a, b, N)

s = Memory{T}(undef, N)
s .= f.(t_range)
#s = Memory{T}([-1.2; 0; 0.5; 1; 1.2; 2; 1;]) # debug.

#c = ITP.get_coeffs(s, ϵ)

# This is the main exported function for this package.
# It is a faster method based on ITP.query, but simplifies the computation for the border 8 samples.

tmp_r = LinRange(a, b, length(s))
#Np = 5 # 3 + 2 # minimum allowed, but might be quite far from samples.
#Np = 8
Np = 10
buf = ITP.FitBuffer1D(T, length(s); N_padding = Np)
#padding_option = ITP.Lagrange4Padding()
#padding_option = ITP.ConstantPadding()
padding_option = ITP.LinearPadding()
#extrapolation_option = ITP.ZeroExtrapolation()
extrapolation_option = ITP.ConstantExtrapolation()
itp1D = ITP.Interpolator1D(
    padding_option,
    extrapolation_option,
    buf, s, a, b; ϵ = ϵ,
) # allocates itp1D.coeffs.

N_params = ITP.get_num_coeffs(itp1D)
#c = randn(rng, T, N_params)
c = f.(LinRange(first(t_range), last(t_range), N_params))
ITP.update_coeffs!(itp1D, c)


# Query range.
query_lb, query_ub = ITP.get_itp_interval(itp1D)

tq_range = LinRange(query_lb, query_ub, Nq)
#tq_range = LinRange(a, b, Nq) # debug for error plot
#tq_range = t_range # debug for error plot

# query.
#query_reference_t = collect( ITP.query(u, c, A) for u in tq_range )
query1D_t = collect(ITP.query1D(u, itp1D) for u in tq_range)

f_tq = f.(tq_range) # oracle

# You can see the border region is not in agreement.
PLT.figure(fig_num)
fig_num += 1
#PLT.plot(t_range, s, label = "data", "x")
#PLT.plot(tq_range, f_tq, label = "oracle")
PLT.plot(tq_range, query1D_t, label = "query1D()")
PLT.legend()

# Test derivative of query with respect to the coefficients.

"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where {T <: Real}
    return (x - a) * (d - c) / (b - a) + c
end

test_lb = a - T(4)
test_ub = b + T(4)
N_tests = 1000
xs = [ convertcompactdomain(rand(rng, T), zero(T), one(T), test_lb, test_ub) for _ in 1:N_tests ]

out = [ ITP.query1D_parameter_derivatives(x, itp1D) for x in xs ]

# q_x, k1, k2, k3, k4, d1, d2, d3, d4

function query1D_params_as_input!(itp, c::AbstractVector, x0::AbstractFloat)
    ITP.update_coeffs!(itp, c)
    return ITP.query1D(x0, itp)
end

zero_tol = T(1.0e-10)
x_test = xs[1]

function test_param_derivatives!(itp1D, c_test, x_test, zero_tol)

    ITP.update_coeffs!(itp1D, c_test)
    h = cc -> query1D_params_as_input!(itp1D, cc, x_test)
    dq_dc_ND = FiniteDiff.finite_difference_gradient(h, c_test)

    q_x_oracle = ITP.query1D(x_test, itp1D)
    q_x, k1, k2, k3, k4, d1, d2, d3, d4 = ITP.query1D_parameter_derivatives(x_test, itp1D)
    # k1 to k4 are the gradient indices that are non-zero.
    # The d1 to d4 values are the gradient values to the corresponding k.

    r = copy(dq_dc_ND)
    r[k1] -= d1
    r[k2] -= d2
    r[k3] -= d3
    r[k4] -= d4
    #@show norm(r) # should be zero if analytic gradient works.
    @assert norm(r) < zero_tol
    return @assert abs(q_x_oracle - q_x) < zero_tol
end

# The following should throw an error if the analytical and numerical derivatives don't match.
N_coeffs_to_test = 50
for _ in 1:N_coeffs_to_test

    randn!(rng, c)

    for x in xs
        test_param_derivatives!(itp1D, c, x, zero_tol)
    end
end

nothing
