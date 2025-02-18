# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

using Test, LinearAlgebra, Random

import CubicBSplineInterpolation as ITP

@testset "1D test" begin

    # test parameters
    T = Float64
    Nq_max = 1213
    zero_tol = eps(T) * 2

    # Set up sample problem.
    ϵ = T(1e-8) # coefficient tolerance.
    N = 100
    a = T(-10)
    b = T(3.4)
    t_range = LinRange(a, b, N)

    f = xx -> exp(-(1 / 17) * (xx / 3)^2) # oracle function.
    s = f.(t_range) # generate interpolation samples.

    N_padding = 5 # This needs to be >= 2 for the interpolation result to fully match up at all sampling locations.
    buf = ITP.FitBuffer1D(T, length(s); N_padding=N_padding)

    # packaged version. Check implementation of interior interval.
    itp1D = ITP.Interpolator1D(buf, s, a, b; ϵ=ϵ)
    lb, ub = ITP.get_itp_interval(itp1D)
    @assert abs(a - lb) < zero_tol
    @assert abs(b - ub) < zero_tol

    # residual fit error should be zero.
    # residual error, should be near zero if N_padding >= 2.
    s_rec = collect(ITP.query1D(u, itp1D) for u in t_range)
    @test norm(s - s_rec) / norm(s) < eps(T)

    # copy-test
    itp3 = deepcopy(itp1D)
    itp2 = ITP.Interpolator1D(itp1D) # deepcopy(.) but should be faster.
    s2 = collect(ITP.query1D(u, itp2) for u in t_range)
    s3 = collect(ITP.query1D(u, itp3) for u in t_range)
    @test norm(s_rec - s2) / norm(s_rec) < eps(T)
    @test norm(s_rec - s3) / norm(s_rec) < eps(T)

    rng = Random.Xoshiro(0)

    # Test LinRange sequences up to length `N_Nq_max`.
    for Nq = N:Nq_max

        query_lb, query_ub = ITP.get_itp_interval(itp1D)
        extrapolation_len = abs(randn(rng, T))
        tq_range = LinRange(
            query_lb - extrapolation_len,
            query_ub + extrapolation_len,
            Nq,
        )
        results = collect(ITP.query1D(u, itp1D) for u in tq_range)

        # sanity check.
        ITP.update_itp!(itp1D, buf, s; ϵ=ϵ) # fits and overwrites existing coefficients.
        results2 = collect(ITP.query1D(u, itp1D) for u in tq_range)
        @test norm(results - results2) / norm(results) < zero_tol # should be zero.
    end
end

import FiniteDiff

"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where {T<:Real}
    return (x - a) * (d - c) / (b - a) + c
end

function query1D_params_as_input!(itp, c::AbstractVector, x0::AbstractFloat)
    ITP.update_coeffs!(itp, c)
    return ITP.query1D(x0, itp)
end

function test_param_derivatives!(itp1D, c_test, x_test, zero_tol)

    ITP.update_coeffs!(itp1D, c_test)
    h = cc -> query1D_params_as_input!(itp1D, cc, x_test)
    dq_dc_ND = FiniteDiff.finite_difference_gradient(h, c_test)

    q_x_oracle = ITP.query1D(x_test, itp1D)
    q_x, k1, k2, k3, k4, d1, d2, d3, d4 = ITP.query1D_parameter_derivatives(x_test, itp1D)

    r = copy(dq_dc_ND)
    r[k1] -= d1
    r[k2] -= d2
    r[k3] -= d3
    r[k4] -= d4
    #@show norm(r) # should be zero if analytic gradient works.
    @assert norm(r) < zero_tol
    @assert abs(q_x_oracle - q_x) < zero_tol
end

# ∂q/∂c, where q is the query function.
@testset "1D coeffs derivative test" begin

    N_coeffs_to_test = 50
    T = Float64
    rng = Random.Xoshiro(0)

    zero_tol = T(1e-10)

    # set up interpolator domain.
    ϵ = T(1e-8)
    N = 1000
    a = T(-10)
    b = T(3.4)
    buf = ITP.FitBuffer1D(T, N; N_padding=10)
    itp1D = ITP.Interpolator1D(
        ITP.LinearPadding(),
        ITP.ConstantExtrapolation(),
        buf, randn(rng, T, N), a, b; ϵ=ϵ,
    )

    test_lb = a - T(4)
    test_ub = b + T(4)
    N_tests = 1000
    xs = [convertcompactdomain(rand(rng, T), zero(T), one(T), test_lb, test_ub) for _ = 1:N_tests]

    N_coeffs = ITP.get_num_coeffs(itp1D)
    c = zeros(T, N_coeffs)
    for _ = 1:N_coeffs_to_test

        randn!(rng, c)

        for x in xs
            test_param_derivatives!(itp1D, c, x, zero_tol)
        end
    end


end

# TODO tests for 2D, 1D complex, 2D complex,