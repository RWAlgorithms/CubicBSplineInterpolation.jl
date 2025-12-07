# SPDX-License-Identifier: MPL-2.0
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

#const T = Float64
T = Float32

ϵ = eps(T) * 2
a1 = T(-4.0)
b1 = T(3.45)
N1 = 100
a2 = T(-2.0)
b2 = T(1.23)
N2 = 113
t_range1 = LinRange(a1, b1, N1)
t_range2 = LinRange(a2, b2, N2)
f = (xx, yy) -> sinc((xx)^2 + (yy)^2 + 20)
S = [f(x1, x2) for x1 in t_range1, x2 in t_range2]


# fit.
buf = ITP.FitBuffer2D(T, size(S); N_padding = (10, 10))
#itp2D = ITP.Interpolator2D(buf, S, a1, b1, a2, b2; ϵ = ϵ)

# padding_option = ITP.LinearPadding()
# extrapolation_option = ITP.ZeroExtrapolation()

padding_option = ITP.LinearPadding()
extrapolation_option = ITP.ConstantExtrapolation()

itp2D = ITP.Interpolator2D(
    padding_option,
    extrapolation_option,
    buf, S, a1, b1, a2, b2; ϵ = ϵ,
)

# test update_itp!
c_back = copy(itp2D.coeffs)
S_random = randn(Random.Xoshiro(0), T, size(S))
ITP.update_itp!(padding_option, extrapolation_option, itp2D, buf, S_random; ϵ = ϵ)
@assert norm(c_back - itp2D.coeffs) > eps(T) * 10

ITP.update_itp!(padding_option, extrapolation_option, itp2D, buf, S; ϵ = ϵ)
@assert norm(c_back - itp2D.coeffs) < eps(T)

# timing
#@btime ITP.update_itp!($itp2D, $buf, $S; ϵ = $ϵ)
"""
size(S) = (100, 113)
250.779 μs (0 allocations: 0 bytes)
"""

#@btime ITP.Interpolator2D($buf, $S, $a1, $b1, $a2, $b2; ϵ = $ϵ)
"""
size(S) = (100, 113)
246.662 μs (4 allocations: 62.50 KiB)
"""

x1 = (t_range1[44] + t_range1[44]) / 2
x2 = (t_range2[44] + t_range2[44]) / 2
out = ITP.query2D(x1, x2, itp2D)
#@btime ITP.query2D($x1, $x2, $itp2D)
"""
32.384 ns (0 allocations: 0 bytes)
"""

# fit residual.
q_S = [ITP.query2D(x1, x2, itp2D) for x1 in t_range1, x2 in t_range2]
println("Residual fit error. This should be near zero.")
@show norm(S - q_S)
@show norm(S - q_S) / norm(S)
println()

# query.

Nq1 = N1 * 10
Nq2 = N2 * 10
aq1, bq1 = t_range1[begin], t_range1[end]
aq2, bq2 = t_range2[begin], t_range2[end]
tq_range1 = LinRange(aq1, bq1, Nq1)
tq_range2 = LinRange(aq2, bq2, Nq2)

Yq = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]
Sq = [f(x1, x2) for x1 in tq_range1, x2 in tq_range2]

println("Real-valued case.")

println("Relative error in the query regions:")
@show norm(Sq - Yq) / norm(Sq)


PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq)
PLT.title("Yq")

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Sq)
PLT.title("Oracle")

# extrapolation
extrapolation_len = T(0.3)
tq_range1 = LinRange(aq1 - extrapolation_len, bq1 + extrapolation_len, Nq1)
tq_range2 = LinRange(aq2 - extrapolation_len, bq2 + extrapolation_len, Nq2)

Yq = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq)
PLT.title("extrapolation length: $(extrapolation_len), Yq")

# @assert 3==4

# Timing
x1 = tq_range1[342]
x2 = tq_range2[432]
#@btime ITP.query2D($x1, $x2, $itp2D);
"""
31.971 ns (0 allocations: 0 bytes)
"""

# Interpolations.jl
import Interpolations

function setup_itp(
        A::Matrix{T},
        A_r,
        A_λ,
    ) where {T <: AbstractFloat}

    real_itp = Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_λ)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    return real_setp
end

itp = setup_itp(S, t_range1, t_range2)

out_itp = itp(x1, x2)
out_query2D = ITP.query2D(x1, x2, itp2D)
out_oracle = f(x1, x2)

@show abs(out_query2D - out_oracle) / norm(out_oracle)
@show abs(out_itp - out_oracle) / norm(out_oracle)

# @btime $f($x1, $x2)
# @btime $itp($x1, $x2) # From Interpolations.jl # does a little worse for Float32 on Ryzen 7 1700.
# @btime ITP.query2D($x1, $x2, $itp2D) # The default method of this package. does a little worse for Float64 on Ryzen 7 1700.
# println()
"""
13.020 ns (0 allocations: 0 bytes)
35.363 ns (0 allocations: 0 bytes)
31.971 ns (0 allocations: 0 bytes)
"""


# # Complex values, 2D

println("Complex-valued case.")

# Specify oracles for the real and imaginary parts.
f_real = (xx, yy) -> sinc((xx)^2 + (yy)^2)
f_imag = (xx, yy) -> tanh((xx)^2 + (yy)^2)
f = (xx, yy) -> Complex(sinc((xx)^2 + (yy)^2), tanh((xx)^2 + (yy)^2)) # This allocates for some reason. Perhaps due to "boxing". Complex(f_real(xx,yy), f_imag(xx,yy))

t_range1 = LinRange(T(-3), T(3), 1000)
t_range2 = LinRange(T(-1.23), T(4.56), 783)

# Generate samples.
Sr = [f_real(x1, x2) for x1 in t_range1, x2 in t_range2]
Si = [f_imag(x1, x2) for x1 in t_range1, x2 in t_range2]

cbuf = ITP.FitBuffer2D(T, size(Sr); N_padding = (8, 7))
citp = ITP.Interpolator2DComplex(cbuf, Sr, Si, first(t_range1), last(t_range1), first(t_range2), last(t_range2))

# test update_itp!
cr_back = copy(citp.real_coeffs)
ci_back = copy(citp.imag_coeffs)
Sr_random = randn(Random.Xoshiro(0), T, size(Sr))
Si_random = randn(Random.Xoshiro(0), T, size(Si))
ITP.update_itp!(citp, cbuf, Sr_random, Si_random; ϵ = ϵ)
@assert norm(cr_back - citp.real_coeffs) > eps(T) * 10
@assert norm(ci_back - citp.imag_coeffs) > eps(T) * 10

ITP.update_itp!(citp, cbuf, Sr, Si; ϵ = ϵ)
@assert norm(cr_back - citp.real_coeffs) < eps(T)
@assert norm(ci_back - citp.imag_coeffs) < eps(T)


# cetermine query intervals.
query_lb1, query_ub1, query_lb2, query_ub2 = ITP.get_itp_interval(citp)
tq_range1 = LinRange(query_lb1, query_ub1, 10000)
tq_range2 = LinRange(query_lb2, query_ub2, 1473)

Yq = [ ITP.query2D(x1, x2, citp) for x1 in tq_range1, x2 in tq_range2 ]
Sq = [ f(x1, x2) for x1 in tq_range1, x2 in tq_range2 ]
@show norm(Sq - Yq) / norm(Sq)

# Compare with Interpolations

function setup_itp(
        S_real::Matrix{T},
        S_imag::Matrix{T},
        A_r,
        A_λ,
    ) where {T <: AbstractFloat}

    real_itp = Interpolations.interpolate(S_real, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_λ)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    imag_itp = Interpolations.interpolate(S_imag, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_λ)
    imag_setp = Interpolations.extrapolate(imag_sitp, zero(T)) # zero outside interp range.

    return real_setp, imag_setp
end

struct ComplexItp{IT}
    real_itp::IT
    imag_itp::IT
end

function query_itp(x, y, A::ComplexItp)
    return Complex(A.real_itp(x, y), A.imag_itp(x, y))
end

itp_real, itp_imag = setup_itp(Sr, Si, t_range1, t_range2)
itp_struct = ComplexItp(itp_real, itp_imag)
itp = (xx, yy) -> query_itp(xx, yy, itp_struct) # doing Complex(itp_real(xx,yy), itp_imag(xx,yy)) leads to allocation for some reason.

out_itp = itp(x1, x2)
out_query2D = ITP.query2D(x1, x2, citp)
out_oracle = f(x1, x2)

@show abs(out_query2D - out_oracle)
@show abs(out_itp - out_oracle)

println("Evaluation timing orders: Oracle, Interpolations.jl, callable wrapper to Interpolations.jl, ITP.query1D")
@btime $f($x1, $x2)
@btime $itp($x1, $x2) # allocates for some reason.
@btime $query_itp($x1, $x2, $itp_struct)
@btime ITP.query2D($x1, $x2, $citp) # The default method of this package. does a little worse for Float64 on Ryzen 7 1700.
println()

versioninfo()

"""
julia> include("demo_2D.jl")
Real-valued case.
Relative error in the interior query regions:
norm(Sq - Yq) / norm(Sq) = 0.000675699783428046
35.413 ns (0 allocations: 0 bytes)
abs(out_query2D - out_oracle) = 4.213215210605026e-5
abs(out_itp - out_oracle) = 4.213215210610577e-5
14.526 ns (0 allocations: 0 bytes)
31.205 ns (0 allocations: 0 bytes)
35.413 ns (0 allocations: 0 bytes)

Complex-valued case.
norm(Sq - Yq) / norm(Sq) = 2.0110752556948815e-6
abs(out_query2D - out_oracle) = 1.1237840788581748e-9
abs(out_itp - out_oracle) = 1.1237841355405732e-9
35.716 ns (0 allocations: 0 bytes)
104.817 ns (3 allocations: 64 bytes)
63.967 ns (0 allocations: 0 bytes)
55.764 ns (0 allocations: 0 bytes)

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
