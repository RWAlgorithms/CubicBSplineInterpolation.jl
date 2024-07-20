
if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

const T = Float64

ϵ = T(1e-16)
Nq = 757
N = 100
a = T(-10)
b = T(3.4)
t = LinRange(a, b, N)
f = sin # the oracle function we use for this test.

A = ITP.IntervalConversion(a, b, N)

s = Memory{T}(undef, N)
s .= f.(t)

c = ITP.get_coeffs(s, ϵ)

#q_t = collect( query(u, c) for u in t )
tq = LinRange(a, b, Nq)
q_t = collect( ITP.query(u, c, A) for u in tq )
q3_t = collect( ITP.query_interior(u, c, A) for u in tq )

# You can use the `Interpolator1D` for a convinence constructor, or construct the query cache `A` and coefficients `c` yourself.
itp1D = ITP.Interpolator1D(s, a, b; ϵ = ϵ)
q4_t = collect( ITP.query_interior(u, itp1D) for u in tq )

@assert norm(q3_t - q4_t) < 1e-14 # both approaches creates the same spline interpolator.


# You can see the border region is not in agreement.
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t, s, label = "data", "x")
PLT.plot(tq, q_t, label = "spline oracle")
PLT.plot(tq, q3_t, label = "spline 3", "--")
PLT.legend()

println("This is the region that does not suffer from border artefacts. The border region is within 4 samples from the edge.")
@show ITP.get_query_interval(itp1D)
println("Use this information to decide your own constraints on whether or not to query the interpolator in the non-border interrior, border, or extrapolation regions.")
println()

query_lb, query_ub = ITP.get_query_interval(itp1D)
println("Relative l-2 error for the the query region, excluding the border 4 samples.")
tq2 = filter(xx->(query_lb < xx < query_ub), tq)
f_tq_exclude_4 = f.(tq2)
q4_t_exclude_4 = collect( ITP.query_interior(u, itp1D) for u in tq2 )
@show norm(f_tq_exclude_4 - q4_t_exclude_4)/norm(f_tq_exclude_4)


# # Compare with Interpolations.jl

_, problem_ind = findmax(abs.(f_tq_exclude_4 - q4_t_exclude_4))
u = tq2[problem_ind]

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

itp = setup_itp(s, t)

println("Discrepancy between Interpolations.jl and our result at our worse fit location of $(u):")
@show abs(itp(u) - ITP.query_interior(u, itp1D))

println("Timing:")
@btime ITP.query_interior($u, $itp1D)
@btime $itp($u)

"""
This is the region that does not suffer from border artefacts. The border region is within 4 samples from the edge.
ITP.get_query_interval(itp1D) = (-9.593939393939394, 2.9939393939393946)
Use this information to decide your own constraints on whether or not to query the interpolator in the non-border interrior, border, or extrapolation regions.

Relative l-2 error for the the query region, excluding the border 4 samples.
norm(f_tq_exclude_4 - q4_t_exclude_4) / norm(f_tq_exclude_4) = 0.0003555412437610131
Discrepancy between Interpolations.jl and our result at our worse fit location of 2.9391534391534386:
abs(itp(u) - ITP.query_interior(u, itp1D)) = 0.003318465340532706
Timing:
    18.107 ns (0 allocations: 0 bytes)
    17.106 ns (0 allocations: 0 bytes)

"""

nothing