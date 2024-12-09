# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>


if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

const T = Float64
#const T = Float32

#ϵ = T(1e-16)
ϵ = T(1e-8)
Nq = 757
N = 1000
#N = 7 # debug
a = T(-10)
b = T(3.4)
t_range = LinRange(a, b, N)
#f = xx->sinc(sin((xx/5)^3+(xx/20)^2+4.32)) # not bad
#f = xx->sin(xx) # not bad.
#f = xx->(exp(-(1/17)*(xx-3)^2) +55) # runge effect near right boundary.
f = xx->exp(-(1/17)*(xx-3)^2) # runge effect near right boundary.
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

## timings
# c_backup = copy(itp1D.coeffs)
# s_SA = view(s, 1:length(s))
# @btime ITP.update_itp!($itp1D, $buf, $s_SA; ϵ = $ϵ); # doesn't allocate.
# @show norm(itp1D.coeffs - c_backup)
# @assert 1==23
"""
N = 1000
5.617 μs (0 allocations: 0 bytes)
"""

#@btime ITP.Interpolator1D($buf, $s, $a, $b; ϵ = $ϵ);
"""
6.115 μs (2 allocations: 8.03 KiB)
"""

# test update_itp!
c_back = copy(itp1D.coeffs)
s_random = randn(Random.Xoshiro(0), T, size(s))
ITP.update_itp!(padding_option, extrapolation_option, itp1D, buf, s_random; ϵ = ϵ)
@assert norm(c_back - itp1D.coeffs) > eps(T)*10

ITP.update_itp!(padding_option, extrapolation_option, itp1D, buf, s; ϵ = ϵ)
@assert norm(c_back - itp1D.coeffs) < eps(T)


# residual error. Should be close to zero.
ts = tmp_r
q_ts = collect( ITP.query1D(u, itp1D) for u in ts )
@show norm(s - q_ts)/norm(s) # 


# Query range.
query_lb, query_ub = ITP.get_itp_interval(itp1D)

tq_range = LinRange(query_lb, query_ub, Nq)
#tq_range = LinRange(a, b, Nq) # debug for error plot
#tq_range = t_range # debug for error plot

# query.
#query_reference_t = collect( ITP.query(u, c, A) for u in tq_range )
query1D_t = collect( ITP.query1D(u, itp1D) for u in tq_range )

f_tq = f.(tq_range) # oracle

# You can see the border region is not in agreement.
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "x")
PLT.plot(tq_range, f_tq, label = "oracle")
PLT.plot(tq_range, query1D_t, label = "query1D()")
PLT.legend()

#@assert 1==23


println("Relative l-2 error over the interpolation region.")
@show norm(f_tq - query1D_t)/norm(f_tq)


# # Test derivative.

f = xx->ITP.query1D(only(xx), itp1D)
df = xx->ITP.query1D_derivative1(xx,itp1D)

import DifferentiationInterface as DI
import Mooncake as MN

rng = Random.Xoshiro(0)
x0 = randn(rng, T, 1)
backend = DI.AutoMooncake(; config=nothing)
prep = DI.prepare_gradient(f, backend, x0)
df_AD = xx->DI.gradient(f, prep, backend, [xx;])[begin]

# ts = tq_range
#ts = LinRange(a-T(3), b+T(2), N*10)
ts = LinRange(T(3.2), T(3.6), N*10)

df_ts_AN = df.(ts)
df_ts_AD = df_AD.(ts)

println("AD: Relative l-2 error over the interpolation region.")
@show norm(df_ts_AN - df_ts_AD)/norm(df_ts_AN)

PLT.figure(fig_num)
fig_num += 1
# PLT.plot(t_range, s, label = "data", "x")
PLT.plot(ts, df_ts_AN, label = "analytical")
PLT.plot(ts, df_ts_AD, label = "AD", "--")
PLT.legend()
PLT.title("Algorithmic differentiation: first derivative of query1D()")

nothing