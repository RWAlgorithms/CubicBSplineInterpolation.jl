
PLT.close("all")
fig_num = 1

const T = Float64

ϵ = T(1e-16)
Nq = 757
N = 100
a = T(-10)
b = T(3.4)
t = LinRange(a, b, N)

A = ITP.IntervalConversion(a, b, N)

s = Memory{T}(undef, N)
s .= sin.(t)

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

# This is the region that does not suffer from border artefacts. The border region is within 4 samples from the edge.
@show ITP.get_query_interval(itp1D)
# Use this information to decide your own constraints on whether or not to query the interpolator in the non-border interrior, border, or extrapolation regions.

nothing