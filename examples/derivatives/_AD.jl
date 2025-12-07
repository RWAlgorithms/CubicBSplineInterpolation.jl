# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# AD via Mooncake.jl is disabled for now.

if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

const T = Float64
#const T = Float32

#ϵ = T(1e-16)
ϵ = T(1.0e-8)
Nq = 17570
N = 1000
#N = 7 # debug
a = T(-10)
b = T(3.4)
t_range = LinRange(a, b, N)
#f = xx->sinc(sin((xx/5)^3+(xx/20)^2+4.32)) # not bad
#f = xx->sin(xx) # not bad.
#f = xx->(exp(-(1/17)*(xx-3)^2) +55)
f = xx -> (sinc(xx) + exp(-(1 / 17) * (xx - 3)^2) + 1)
f = xx -> (exp(-(1 / 17) * (xx - 3)^2) + 1)
f = xx -> (exp(-(1 / 17) * (xx - 3)^2) + sin(xx - 1.2))

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
padding_option = ITP.Lagrange4Padding()
#padding_option = ITP.ConstantPadding()
#padding_option = ITP.LinearPadding()
#extrapolation_option = ITP.ZeroExtrapolation()
extrapolation_option = ITP.ConstantExtrapolation()
itp1D = ITP.Interpolator1D(
    padding_option,
    extrapolation_option,
    buf, s, a, b; ϵ = ϵ,
) # allocates itp1D.coeffs.


# test update_itp!
c_back = copy(itp1D.coeffs)
s_random = randn(Random.Xoshiro(0), T, size(s))
ITP.update_itp!(padding_option, extrapolation_option, itp1D, buf, s_random; ϵ = ϵ)
@assert norm(c_back - itp1D.coeffs) > eps(T) * 10

ITP.update_itp!(padding_option, extrapolation_option, itp1D, buf, s; ϵ = ϵ)
@assert norm(c_back - itp1D.coeffs) < eps(T)


# residual error. Should be close to zero.
ts = tmp_r
q_ts = collect(ITP.query1D(u, itp1D) for u in ts)
@show norm(s - q_ts) / norm(s) #


# Query range.
query_lb, query_ub = ITP.get_itp_interval(itp1D)

#tq_range = LinRange(query_lb, query_ub, Nq)
tq_range = LinRange(query_lb - 1, query_ub + 1, Nq)
tq_range = LinRange(query_lb - 10, query_ub + 10, Nq)

## zoom.
#tq_range = LinRange(T(query_ub-0.1), T(query_ub+0.1), Nq)
#tq_range = LinRange(T(3.3), T(3.5), Nq)

# query.
query1D_t = collect(ITP.query1D(u, itp1D) for u in tq_range)
#query1DLE_t = collect( ITP.query1D_linear_extrapolation(u, itp1D) for u in tq_range )


etp1D = ITP.QuadraticExtrapolator1D(itp1D)
q_tq = collect(ITP.query1D(u, etp1D) for u in tq_range)


f_tq = f.(tq_range) # oracle

# You can see the border region is not in agreement.
PLT.figure(fig_num)
fig_num += 1
#PLT.plot(t_range, s, label = "data", "x")
PLT.plot(tq_range, f_tq, label = "oracle")
PLT.plot(tq_range, query1D_t, label = "query1D()", "--")
#PLT.plot(tq_range, query1DLE_t, label = "linear extrapolation", "--")
PLT.plot(tq_range, q_tq, label = "extrapolation", "x")
PLT.legend()
PLT.title("Extrapolation variant of query1D()")

# view derivatives

df_tqs = [ ITP.query1D_derivative1(u, itp1D) for u in tq_range ]
d2f_tqs = [ ITP.query1D_derivative2(u, itp1D) for u in tq_range ]

PLT.figure(fig_num)
fig_num += 1
PLT.plot(tq_range, query1D_t, label = "query function")
PLT.plot(tq_range, df_tqs, label = "first derivative", "--")
PLT.plot(tq_range, d2f_tqs, label = "second derivative", "x")
PLT.legend()
PLT.title("Derivatives")
# Second derivative is not differentiable


println("Relative l-2 error over the interpolation region.")
@show norm(f_tq - query1D_t) / norm(f_tq)


# # Test derivative.

f = xx -> ITP.query1D(only(xx), itp1D)
df = xx -> ITP.query1D_derivative1(xx, itp1D)

# import DifferentiationInterface as DI
# import Mooncake as MN

# rng = Random.Xoshiro(0)
# x0 = randn(rng, T, 1)
# backend = DI.AutoMooncake(; config = nothing)
# prep = DI.prepare_gradient(f, backend, x0)
# df_AD = xx -> DI.gradient(f, prep, backend, [xx;])[begin]

# # ts = tq_range
# #ts = LinRange(a-T(3), b+T(2), N*10)
# ts = LinRange(T(3.2), T(3.6), N * 10)

# df_ts_AN = df.(ts)
# df_ts_AD = df_AD.(ts)

# println("AD: Relative l-2 error over the interpolation region.")
# @show norm(df_ts_AN - df_ts_AD) / norm(df_ts_AN)

# PLT.figure(fig_num)
# fig_num += 1
# # PLT.plot(t_range, s, label = "data", "x")
# PLT.plot(ts, df_ts_AN, label = "analytical")
# PLT.plot(ts, df_ts_AD, label = "AD", "--")
# PLT.legend()
# PLT.title("Algorithmic differentiation: first derivative of query1D()")

nothing
