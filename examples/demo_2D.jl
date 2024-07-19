# # 2D

PLT.close("all")
fig_num = 1

const T = Float64

ϵ = eps(T)*100
a1 = -4.0
b1 = 3.45
N1 = 100
a2 = -2.0
b2 = 1.23
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

#query_interior(xx, yy, C, A1, A2)
Yq = [ITP.query_interior(x1, x2, C, A1, A2) for x1 in tq_range1, x2 in tq_range2]
Sq = [f(x1, x2) for x1 in tq_range1, x2 in tq_range2]

println("Relative error in the interior query regions:")
@show norm(Sq-Yq)/norm(Sq)

# test packaged version.
itp2D = ITP.Interpolator2D(S, a1, b1, a2, b2; ϵ = ϵ)
Yq_pkg = [ITP.query_interior(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]
@assert norm(Yq_pkg - Yq) < 1e-14


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
@btime ITP.query_interior($x1, $x2, $itp2D);


# Interpolations.jl
import Interpolations

function setupclpartitionitp(
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

itp = setupclpartitionitp(S, t_range1, t_range2)

out_itp = itp(x1, x2)
out_mine = ITP.query_interior(x1, x2, itp2D)
out_oracle = f(x1,x2)

@show abs(out_mine - out_oracle)
@show abs(out_itp - out_oracle)

@btime $f($x1, $x2)
@btime $itp($x1, $x2)
@btime ITP.query_interior($x1, $x2, $itp2D)

"""
julia> include("demo_2D.jl")
Relative error in the interior query regions:
norm(Sq - Yq) / norm(Sq) = 0.0006756997834280465
43.253 ns (0 allocations: 0 bytes)
abs(out_mine - out_oracle) = 4.213215210605026e-5
abs(out_itp - out_oracle) = 4.213215210610577e-5
14.506 ns (0 allocations: 0 bytes)
30.822 ns (0 allocations: 0 bytes)
42.948 ns (0 allocations: 0 bytes)

# with N1 and N2, x10 each.
julia> include("demo_2D.jl")
Relative error in the interior query regions:
norm(Sq - Yq) / norm(Sq) = 8.299799743547375e-5
43.253 ns (0 allocations: 0 bytes)
abs(out_mine - out_oracle) = 9.777807296468266e-9
abs(out_itp - out_oracle) = 9.777807292998819e-9
15.058 ns (0 allocations: 0 bytes)
30.853 ns (0 allocations: 0 bytes)
43.253 ns (0 allocations: 0 bytes)

"""

nothing