
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
N = 100
a = T(-10)
b = T(3.4)
t_range = LinRange(a, b, N)
#f = xx->sinc(sin(xx^3+xx^2+4.32)) # the oracle function we use for this test.
#f = xx->sin(xx)
f = xx->exp(-(1/17)*(xx-3)^2)

A = ITP.IntervalConversion(a, b, N)

s = Memory{T}(undef, N)
s .= f.(t_range)

c = ITP.get_coeffs(s, ϵ)

# This is the main exported function for this package.
# It is a faster method based on ITP.query, but simplifies the computation for the border 8 samples.
itp1D = ITP.Interpolator1D(s, a, b; ϵ = ϵ)

# Query range.
query_lb, query_ub = ITP.get_query_interval(itp1D)
@assert isapprox(t_range[8], query_lb)
@assert isapprox(t_range[end-7], query_ub) # not exact to do floating precision issues.

tq_range = LinRange(query_lb, query_ub, Nq)
#tq_range = LinRange(a, b, Nq) # debug for error plot
#tq_range = t_range # debug for error plot

# query.
query_reference_t = collect( ITP.query(u, c, A) for u in tq_range )
query1D_t = collect( ITP.query1D(u, itp1D) for u in tq_range )


# You can see the border region is not in agreement.
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "x")
PLT.plot(tq_range, f.(tq_range), label = "oracle")
PLT.plot(tq_range, query_reference_t, label = "query()")
PLT.plot(tq_range, query1D_t, label = "query1D()", "--")
PLT.legend()


println("This is the region that does not suffer from border artefacts. The border region is within 8 samples from the edge.")
@show ITP.get_query_interval(itp1D)
println("Use this information to decide your own constraints on whether or not to query the interpolator in the non-border interrior, border, or extrapolation regions.")
println()


println("Relative l-2 error for the the query region, excluding the border 8 samples.")
tq2 = filter(xx->(query_lb < xx < query_ub), tq_range)
f_tq_exclude_4 = f.(tq2)
query1D_t_exclude_4 = collect( ITP.query1D(u, itp1D) for u in tq2 )
@show norm(f_tq_exclude_4 - query1D_t_exclude_4)/norm(f_tq_exclude_4)


# # Compare with Interpolations.jl

_, problem_ind = findmax(abs.(f_tq_exclude_4 - query1D_t_exclude_4))
xq = tq2[problem_ind]

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

println("Timing:")
@btime ITP.query1D($xq, $itp1D)
@btime $itp($xq)

itp_tq = itp.(tq_range)

f_tq = f.(tq_range)
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label = "data", "x")
PLT.plot(tq_range, log.(abs.(query1D_t - f_tq)), label = "query1D: log error")
PLT.plot(tq_range, log.(abs.(itp_tq - f_tq)), "o", label = "Interpolations.jl: log error")
PLT.plot(tq_range, log.(abs.(query_reference_t - f_tq)), "--", label = "query: log error")
PLT.title("Error plot")
PLT.legend()

# There can be a bit of error for the region in the border 8 samples.
coeffs_itp = itp.itp.itp

PLT.figure(fig_num)
fig_num += 1
PLT.plot(abs.(itp1D.coeffs - coeffs_itp), "x")
PLT.title("Coefficients difference between Interpolations.jl and the one in this package.")
PLT.legend()


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
citp = ITP.Interpolator1DComplex(Sr, Si, first(t_range), last(t_range))

# cetermine query intervals.
query_lb, query_ub = ITP.get_query_interval(citp)
tq_range = LinRange(query_lb, query_ub, 10000)

Yq = [ ITP.query1D(x1, citp) for x1 in tq_range ]
Sq = [ f(x1) for x1 in tq_range ]
@show norm(Sq-Yq)/norm(Sq)

# Choose the place where we have the most discrepancy as our single query point.
_, problem_index = findmax(abs.(Sq-Yq))
xq = tq_range[problem_index]

# Compare with Interpolations

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

@btime $f($xq)
@btime $itp($xq) # allocates for some reason.
@btime $query_itp($xq, $itp_struct)
@btime ITP.query1D($xq, $citp) # The default method of this package. does a little worse for Float64 on Ryzen 7 1700.
println()

versioninfo()

"""
julia> include("demo_1D.jl")
This is the region that does not suffer from border artefacts. The border region is within 8 samples from the edge.
ITP.get_query_interval(itp1D) = (-9.052525252525252, 2.4525252525252528)
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