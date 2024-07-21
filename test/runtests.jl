using Test, LinearAlgebra

import CubicBSplineInterpolation as ITP

@testset "1D test" begin

    # test parameters
    T = Float64
    Nq_max = 1213
    zero_tol = eps(T)*100

    # Set up sample problem.
    ϵ = eps(T)*1_000 # coefficient tolerance.
    N = 100
    a = T(-10)
    b = T(3.4)
    t = LinRange(a, b, N)

    A = ITP.IntervalConversion(a, b, N)

    s = Memory{T}(undef, N)
    s .= sin.(t) # oracle function.

    c = ITP.get_coeffs(s, ϵ)
    
    # Valid query range.
    aq = t[4]
    bq = t[end-3]
    
    # packaged version. Check implementation of interior interval.
    itp1D = ITP.Interpolator1D(s, a, b; ϵ = ϵ)
    lb, ub = ITP.get_query_interval(itp1D)
    @assert abs(aq - lb) < zero_tol
    @assert abs(bq - ub) < zero_tol

    # Test LinRange sequences up to length `N_Nq_max`.
    for Nq = N:Nq_max
        tq = LinRange(aq, bq, Nq)
        q_t = collect( ITP.query(u, c, A) for u in tq )
        q3_t = collect( ITP.query1D(u, itp1D) for u in tq )
        
        seq_test = q3_t[begin+3:end-3]
        seq_oracle = q_t[begin+3:end-3]
        @test norm(seq_test-seq_oracle)/norm(seq_oracle) < zero_tol
    end
end