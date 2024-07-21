if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

# # 2D


PLT.close("all")
fig_num = 1

#const T = Float64
const T = Float32

ϵ = eps(T)*100
a1 = T(-4.0)
b1 = T(3.45)
N1 = 100
a2 = T(-2.0)
b2 = T(1.23)
N2 = 113
t_range1 = LinRange(a1, b1, N1)
t_range2 = LinRange(a2, b2, N2)
f = (xx,yy)->sinc((xx)^2+(yy)^2)
S = [f(x1,x2) for x1 in t_range1, x2 in t_range2]
C = ITP.get_coeffs(S, ϵ)

A1 = ITP.IntervalConversion(a1, b1, N1)
A2 = ITP.IntervalConversion(a2, b2, N2)

# query
Nq1 = N1*10
Nq2 = N2*10
aq1, bq1 = t_range1[4], t_range1[end-3]
aq2, bq2 = t_range2[4], t_range2[end-3]
tq_range1 = LinRange(aq1, bq1, Nq1)
tq_range2 = LinRange(aq2, bq2, Nq2)

# Query spline model and oracle.
itp2D = ITP.Interpolator2D(S, a1, b1, a2, b2; ϵ = ϵ)
Yq = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]
Sq = [f(x1, x2) for x1 in tq_range1, x2 in tq_range2]

println("Real-valued case.")

println("Relative error in the interior query regions:")
@show norm(Sq-Yq)/norm(Sq)


PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq)
PLT.title("Yq")

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Sq)
PLT.title("Oracle")

# Timing
x1 = tq_range1[342]
x2 = tq_range2[432]
@btime ITP.query2D($x1, $x2, $itp2D);


# Interpolations.jl
import Interpolations

function setup_itp(
    A::Matrix{T},
    A_r,
    A_λ,
    ) where T <: AbstractFloat

    real_itp = Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_λ)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    return real_setp
end

itp = setup_itp(S, t_range1, t_range2)

out_itp = itp(x1, x2)
out_query2D = ITP.query2D(x1, x2, itp2D)
out_oracle = f(x1,x2)

@show abs(out_query2D - out_oracle)
@show abs(out_itp - out_oracle)

@btime $f($x1, $x2)
@btime $itp($x1, $x2) # From Interpolations.jl # does a little worse for Float32 on Ryzen 7 1700.
@btime ITP.query2D($x1, $x2, $itp2D) # The default method of this package. does a little worse for Float64 on Ryzen 7 1700.
println()


# # Complex values, 2D

println("Complex-valued case.")

# Specify oracles for the real and imaginary parts.
f_real = (xx,yy)->sinc((xx)^2+(yy)^2)
f_imag = (xx,yy)->tanh((xx)^2+(yy)^2)
f = (xx,yy)->Complex(sinc((xx)^2+(yy)^2), tanh((xx)^2+(yy)^2)) # This allocates for some reason. Perhaps due to "boxing". Complex(f_real(xx,yy), f_imag(xx,yy))

t_range1 = LinRange(T(-3), T(3), 1000)
t_range2 = LinRange(T(-1.23), T(4.56), 783)

# Generate samples.
Sr = [f_real(x1,x2) for x1 in t_range1, x2 in t_range2]
Si = [f_imag(x1,x2) for x1 in t_range1, x2 in t_range2]

citp = ITP.Interpolator2DComplex(Sr, Si, first(t_range1), last(t_range1), first(t_range2), last(t_range2))

# cetermine query intervals.
query_lb1, query_ub1, query_lb2, query_ub2 = ITP.get_query_interval(citp)
tq_range1 = LinRange(query_lb1, query_ub1, 10000)
tq_range2 = LinRange(query_lb2, query_ub2, 1473)


Yq = [ ITP.query2D(x1, x2, citp) for x1 in tq_range1, x2 in tq_range2 ]
Sq = [ f(x1, x2) for x1 in tq_range1, x2 in tq_range2 ]
@show norm(Sq-Yq)/norm(Sq)

# Compare with Interpolations

function setup_itp(
    S_real::Matrix{T},
    S_imag::Matrix{T},
    A_r,
    A_λ,
    ) where T <: AbstractFloat

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
    return Complex(A.real_itp(x,y), A.imag_itp(x,y))
end

itp_real, itp_imag = setup_itp(Sr, Si, t_range1, t_range2)
itp_struct = ComplexItp(itp_real, itp_imag)
itp = (xx,yy)->query_itp(xx,yy, itp_struct) # doing Complex(itp_real(xx,yy), itp_imag(xx,yy)) leads to allocation for some reason.

out_itp = itp(x1, x2)
out_query2D = ITP.query2D(x1, x2, citp)
out_oracle = f(x1,x2)

@show abs(out_query2D - out_oracle)
@show abs(out_itp - out_oracle)

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