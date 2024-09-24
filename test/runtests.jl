using Test, LinearAlgebra, Random

import CubicBSplineInterpolation as ITP

@testset "1D test" begin

    # test parameters
    T = Float64
    Nq_max = 1213
    zero_tol = eps(T)*2

    # Set up sample problem.
    ϵ = T(1e-8) # coefficient tolerance.
    N = 100
    a = T(-10)
    b = T(3.4)
    t_range = LinRange(a, b, N)

    f = xx->exp(-(1/17)*(xx/3)^2) # oracle function.
    s = f.(t_range) # generate interpolation samples.

    N_padding = 2 # This needs to be >= 2 for the interpolation result to fully match up at all sampling locations.
    buf = ITP.FitBuffer1D(T, length(s); N_padding = N_padding)
    
    # packaged version. Check implementation of interior interval.
    itp1D = ITP.Interpolator1D(buf, s, a, b; ϵ = ϵ)
    lb, ub = ITP.get_itp_interval(itp1D)
    @assert abs(a - lb) < zero_tol
    @assert abs(b - ub) < zero_tol

    # residual fit error should be zero.
    # residual error, should be near zero if N_padding >= 2.
    s_rec = collect( ITP.query1D(u, itp1D) for u in t_range )
    @test norm(s-s_rec)/norm(s) < eps(T)

    rng = Random.Xoshiro(0)

    # Test LinRange sequences up to length `N_Nq_max`.
    for Nq = N:Nq_max
        
        query_lb, query_ub = ITP.get_itp_interval(itp1D)
        extrapolation_len = abs(randn(rng, T))
        tq_range = LinRange(
            query_lb- extrapolation_len,
            query_ub + extrapolation_len,
            Nq,
        )
        results = collect( ITP.query1D(u, itp1D) for u in tq_range )
        
        # sanity check.
        ITP.update_itp!(itp1D, buf, s; ϵ = ϵ) # fits and overwrites existing coefficients.
        results2 = collect( ITP.query1D(u, itp1D) for u in tq_range )
        @test norm(results - results2)/norm(results) < zero_tol # should be zero.
    end
end