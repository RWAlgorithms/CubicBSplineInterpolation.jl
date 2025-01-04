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

#@assert 2 ==23

# # Compare with Interpolations.jl

_, problem_ind = findmax(abs.(f_tq- query1D_t))
xq = tq_range[problem_ind]

# Interpolations.jl
import Interpolations

function setup_itp(A_in::Memory{T}, A_r) where T <: AbstractFloat

    A = Array(A_in)
    real_itp = Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    return real_setp
end

itp = setup_itp(s, t_range)

println("Discrepancy between Interpolations.jl and our result at our worse fit location of `xq` = $(xq):")
@show abs(itp(xq) - ITP.query1D(xq, itp1D))/abs(itp(xq))

# println("Query timing:")
# @btime ITP.query1D($xq, $itp1D)
# @btime $itp($xq)
# println()

# println("Setup timing:")
# @btime ITP.Interpolator1D($buf, $s, $a, $b; ϵ = $ϵ)

# s_array = Array(s)
# @btime Interpolations.interpolate($s_array, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
# println()


"""
Setup timing:
262.652 μs (0 allocations: 0 bytes)
3.314 ms (78 allocations: 4.96 MiB)
"""

itp_tq = itp.(tq_range)

f_tq = f.(tq_range)


PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "o")
PLT.plot(tq_range, f_tq, label = "oracle")
PLT.plot(tq_range, query1D_t, label = "query1D")
PLT.plot(tq_range, itp_tq, "--", label = "Interpolations.jl")
PLT.title("Interpolation results")
PLT.legend()

# TODO I am here. look into why Interpolation's coeffs are rock solid.

PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "x")
PLT.plot(tq_range, log.(abs.(query1D_t - f_tq)), label = "query1D: log error")
PLT.plot(tq_range, log.(abs.(itp_tq - f_tq)), "o", label = "Interpolations.jl: log error")
#PLT.plot(tq_range, log.(abs.(query_reference_t - f_tq)), "--", label = "query: log error")
PLT.title("Error plot")
PLT.legend()

# # There can be a bit of error for the region in the border 8 samples.
# coeffs_itp = itp.itp.itp
# coeff_diff = abs.(itp1D.coeffs[Np+1:end-Np] - coeffs_itp)
# PLT.figure(fig_num)
# fig_num += 1
# PLT.plot(coeff_diff, "x")
# PLT.title("Coefficients difference between Interpolations.jl and the one in this package.")


# ## Extrapolaion

#extrapolation_len = T(0.2)
extrapolation_len = T(20.0)
Nq = 15000
tq_range = LinRange(query_lb - extrapolation_len, query_ub + extrapolation_len, Nq)


itp_tq = itp.(tq_range)

# M0 = Np-1
# v = view(itp1D.coeffs, 1:M0)
# #fill!(v, 0)
# tmp = itp1D.coeffs[M0]
# fill!(v, tmp)

# v = view(itp1D.coeffs, (length(itp1D.coeffs)-M0+1):length(itp1D.coeffs))
# #fill!(v, 0)
# tmp = itp1D.coeffs[(length(itp1D.coeffs)-M0+1)]
# fill!(v, tmp)

query1D_t = collect( ITP.query1D(u, itp1D) for u in tq_range )
f_tq = f.(tq_range)


PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "o")
PLT.plot(tq_range, f_tq, label = "oracle")
PLT.plot(tq_range, query1D_t, label = "query1D")
PLT.plot(tq_range, itp_tq, "--", label = "Interpolations.jl")
PLT.title("Extrapolation results")
PLT.legend()

# First derivative
mq = tq_range
dq_tq = [ ITP.query1D_derivative1(u, itp1D) for u in mq ]

h = xx->ITP.query1D(xx, itp1D)
dq_tq_ND = [ FiniteDiff.finite_difference_derivative(h, u) for u in mq ]

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq, dq_tq, label = "implemented")
PLT.plot(mq, dq_tq_ND, "--", label = "numerical")
PLT.title("First derivatives: numerical vs analytical")
PLT.legend()

println("First derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(dq_tq - dq_tq_ND)/norm(dq_tq_ND)

# Second derivative
mq = tq_range
d2q_tq = [ ITP.query1D_derivative2(u, itp1D) for u in mq ]

dh = xx->FiniteDiff.finite_difference_derivative(h, xx) 
d2q_tq_ND = [ FiniteDiff.finite_difference_derivative(dh, u) for u in mq ]

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq, d2q_tq, label = "implemented")
PLT.plot(mq, d2q_tq_ND, "--", label = "numerical")
PLT.title("Second derivatives: numerical vs analytical")
PLT.legend()

println("Second derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(d2q_tq - d2q_tq_ND)/norm(d2q_tq_ND)
println()

# # Complex values, 1D

println("Complex-valued case.")

# Specify oracles for the real and imaginary parts.
f_real = (xx)->sinc((xx)^2)
f_imag = (xx)->tanh((xx)^2)
f = (xx)->Complex(sinc((xx)^2), tanh((xx)^2)) # This allocates for some reason. Perhaps due to "boxing". Complex(f_real(xx,yy), f_imag(xx,yy))

t_range = LinRange(T(-3), T(3), 1000)

# Generate samples.
Sr = [f_real(x1) for x1 in t_range]
Si = [f_imag(x1) for x1 in t_range]

# create interpolator model.
cNp = 5 # the minimum allowed.
#cNp = 10 # smoother boundary.
cbuf = ITP.FitBuffer1D(T, length(Sr); N_padding = cNp)
citp = ITP.Interpolator1DComplex(cbuf, Sr, Si, first(t_range), last(t_range); ϵ = ϵ)


# test update_itp!
cr_back = copy(citp.real_coeffs)
ci_back = copy(citp.imag_coeffs)
Sr_random = randn(Random.Xoshiro(0), T, size(Sr))
Si_random = randn(Random.Xoshiro(0), T, size(Si))
ITP.update_itp!(citp, cbuf, Sr_random, Si_random; ϵ = ϵ)
@assert norm(cr_back - citp.real_coeffs) > eps(T)*10
@assert norm(ci_back - citp.imag_coeffs) > eps(T)*10

ITP.update_itp!(citp, cbuf, Sr, Si; ϵ = ϵ)
@assert norm(cr_back - citp.real_coeffs) < eps(T)
@assert norm(ci_back - citp.imag_coeffs) < eps(T)


# check interpolation fit.
qc = xx->ITP.query1D(xx, citp)
qc_s = qc.(t_range)
println("Interpolation fit residual:")
@show norm(Sr .+ im .* Si - qc_s) # TODO make into test.
println()


# determine query intervals.
query_lb, query_ub = ITP.get_itp_interval(citp)
tq_range = LinRange(query_lb-0.1, query_ub+0.1, 10000)

Yq = [ ITP.query1D(x1, citp) for x1 in tq_range ]
Sq = [ f(x1) for x1 in tq_range ]
@show norm(Sq-Yq)/norm(Sq)

# Choose the place where we have the most discrepancy as our single query point.
_, problem_index = findmax(abs.(Sq-Yq))
xq = tq_range[problem_index]

# visualize.
mq = tq_range
q_tq_real = [ real(ITP.query1D(u, citp)) for u in mq ]
q_tq_imag = [ imag(ITP.query1D(u, citp)) for u in mq ]

PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, Sr, "o", label = "real, data")
PLT.plot(t_range, Si, "x", label = "imag, data")
PLT.plot(mq, q_tq_real, label = "implemented, real")
PLT.plot(mq, q_tq_imag, label = "implemented, imag")
PLT.legend()
PLT.title("interpolationr esult for complex-valued case.")
#@assert 45==5
# First derivative
println("Complex-valued case:")
mq = tq_range
dq_tq_real = [ ITP.query1D_derivative1(u, citp)[1] for u in mq ]
dq_tq_imag = [ ITP.query1D_derivative1(u, citp)[2] for u in mq ]

hr = xx->real(ITP.query1D(xx, citp))
hi = xx->imag(ITP.query1D(xx, citp))
dq_tq_ND_real = [ FiniteDiff.finite_difference_derivative(hr, u) for u in mq ]
dq_tq_ND_imag = [ FiniteDiff.finite_difference_derivative(hi, u) for u in mq ]

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq, dq_tq_real, label = "implemented, real")
PLT.plot(mq, dq_tq_imag, label = "implemented, imag")
PLT.plot(mq, dq_tq_ND_real, "o", label = "numerical, real")
PLT.plot(mq, dq_tq_ND_imag, "x", label = "numerical, imag")
PLT.title("First derivatives: numerical vs analytical")
PLT.legend()

println("First derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(dq_tq_real - dq_tq_ND_real)/norm(dq_tq_ND_real)
@show norm(dq_tq_imag - dq_tq_ND_imag)/norm(dq_tq_ND_imag)

# Second derivative
mq = tq_range
d2q_tq_real = [ ITP.query1D_derivative2(u, citp)[1] for u in mq ]
d2q_tq_imag = [ ITP.query1D_derivative2(u, citp)[2] for u in mq ]

# # hr = xx->ITP.query1D_derivative1(xx, citp)[1]
# # hi = xx->ITP.query1D_derivative1(xx, citp)[2]
# dhr = xx->FiniteDiff.finite_difference_derivative(hr, xx) 
# dhi = xx->FiniteDiff.finite_difference_derivative(hi, xx) 
# d2q_tq_ND_real = [ FiniteDiff.finite_difference_derivative(dhr, u) for u in mq ]
# d2q_tq_ND_imag = [ FiniteDiff.finite_difference_derivative(dhi, u) for u in mq ]
d2q_tq_ND_real = [ FiniteDiff.finite_difference_hessian(uu->hr(uu[1]), [u;])[1] for u in mq ]
d2q_tq_ND_imag = [ FiniteDiff.finite_difference_hessian(uu->hi(uu[1]), [u;])[1] for u in mq ]


PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq, d2q_tq_real, label = "implemented, real")
PLT.plot(mq, d2q_tq_imag, label = "implemented, imag")
PLT.plot(mq, d2q_tq_ND_real, "o", label = "numerical, real")
PLT.plot(mq, d2q_tq_ND_imag, "x", label = "numerical, imag")
PLT.title("Second derivatives: numerical vs analytical")
PLT.legend()

println("Second derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(d2q_tq_real - d2q_tq_ND_real)/norm(d2q_tq_ND_real)
@show norm(d2q_tq_imag - d2q_tq_ND_imag)/norm(d2q_tq_ND_imag)
println()

# ## Compare with Interpolations

function setup_itp(
    S_real::Union{Vector{T}, Memory{T}},
    S_imag::Union{Vector{T}, Memory{T}},
    A_r,
    ) where T <: AbstractFloat

    real_itp = Interpolations.interpolate(S_real, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    imag_itp = Interpolations.interpolate(S_imag, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r)
    imag_setp = Interpolations.extrapolate(imag_sitp, zero(T)) # zero outside interp range.

    return real_setp, imag_setp
end

struct ComplexItp{IT}
    real_itp::IT
    imag_itp::IT
end

function query_itp(x, A::ComplexItp)
    return Complex(A.real_itp(x), A.imag_itp(x))
end

itp_real, itp_imag = setup_itp(Sr, Si, t_range)
itp_struct = ComplexItp(itp_real, itp_imag)
itp = (xx)->query_itp(xx, itp_struct) # doing Complex(itp_real(xx,yy), itp_imag(xx,yy)) leads to allocation for some reason.

out_itp = itp(xq)
out_query1D = ITP.query1D(xq, citp)
out_oracle = f(xq)

@show abs(out_query1D - out_oracle)
@show abs(out_itp - out_oracle)


println("Evaluation timing orders: Oracle, Interpolations.jl, callable wrapper to Interpolations.jl, ITP.query1D")
@btime $f($xq)
@btime $itp($xq) # allocates for some reason.
@btime $query_itp($xq, $itp_struct)
@btime ITP.query1D($xq, $citp) # The default method of this package. does a little worse for Float64 on Ryzen 7 1700.
println()

versioninfo()

"""
julia> include("demo_1D.jl")
This is the region that does not suffer from border artefacts. The border region is within 8 samples from the edge.
ITP.get_itp_interval(itp1D) = (-9.052525252525252, 2.4525252525252528)
Use this information to decide your own constraints on whether or not to query the interpolator in the non-border interrior, border, or extrapolation regions.

Relative l-2 error for the the query region, excluding the border 8 samples.
norm(f_tq_exclude_4 - query1D_t_exclude_4) / norm(f_tq_exclude_4) = 1.0769878543338694e-5
Discrepancy between Interpolations.jl and our result at our worse fit location of `xq` = 2.4068702902036243:
abs(itp(xq) - ITP.query1D(xq, itp1D)) / abs(itp(xq)) = 5.838610316292021e-5
Timing:
17.477 ns (0 allocations: 0 bytes)
17.136 ns (0 allocations: 0 bytes)
No artists with labels found to put in legend.  Note that artists whose label start with an underscore are ignored when legend() is called with no argument.
Complex-valued case.
norm(Sq - Yq) / norm(Sq) = 1.4161575373924584e-6
abs(out_query1D - out_oracle) = 5.836378257648921e-5
abs(out_itp - out_oracle) = 9.870535198214837e-9
29.693 ns (0 allocations: 0 bytes)
85.670 ns (2 allocations: 48 bytes)
42.089 ns (0 allocations: 0 bytes)
21.484 ns (0 allocations: 0 bytes)

Julia Version 1.11.0-rc1
Commit 3a35aec36d1 (2024-06-25 10:23 UTC)
Build Info:
Official https://julialang.org/ release
Platform Info:
OS: Linux (x86_64-linux-gnu)
CPU: 16 × AMD Ryzen 7 1700 Eight-Core Processor
WORD_SIZE: 64
LLVM: libLLVM-16.0.6 (ORCJIT, znver1)
Threads: 1 default, 0 interactive, 1 GC (on 16 virtual cores)
"""

nothing