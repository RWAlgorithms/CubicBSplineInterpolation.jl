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
N = 1000 # Number of interpolation samples.

# This is the oracle function for this example.
f = xx -> (cos(xx) * exp(-(1 / 17) * (xx - 3)^2) + 55)

# `t_range` is the time stamp for the interpolation samples.
a = T(-2) # lower boundary.
b = T(3.4) # upper boundary.
t_range = LinRange(a, b, N)

# `s` is the interpolation samples for this example.
s = Memory{T}(f.(t_range))

# # Boundary oscillation control
# `Since the algorithm used in this package was designed to interpolate non-compactly supported discrete signals, Runge-like oscillations can occur when a compactly supported discrete signal, like `s` is used. Instead of implementing various boundary conditions like in most mainstream spline interpolation libraries, this package opts to control oscillations at the boundaries by padding samples to the input sample set `s`. These additional samples are called *padded samples* throughout this package.

# `Np` is the amount of padding added to each boundary, measured in number of samples. For cubic b-spline interpolation, it has to be 5 or greater. A smaller number means less padded samples to compute in the interpolator constructors, but the chance for larger runge-like oscillations near the 
Np = 10
buf = ITP.FitBuffer1D(T, length(s); N_padding=Np)

# `padding_option` controls the presence of Runge-like oscillations near the boundaries. It is a trait struct with supertype `PaddingOption`. They specify the extrapolation algorithm that computes the padded samples. These algorithms use a number of samples that are closest to the boundary. The possible struct constructors are:
# - `Lagrange4Padding()`: Lagrange extrapolation based on the 4 boundary samples.
# - `ConstantPadding`: constant extrapolation with the constant being the boundary sample.
# - `LinearPadding`: linear extrapolation based on the 2 boundary samples.
padding_option = ITP.LinearPadding()

# # Extrapolation Specification
# `extrapolation_option` is a trait struct with supertype `ExtrapolationOption` that specifies the extrapolation behavior. The possible struct constructors are:
# - `ZeroExtrapolation`: constant extrapolation on the boundary input sample from `s`, i.e., not the boundary samples from the padded samples. 
# - `ZeroExtrapolation`: constant extrapolation where the constant is zero.
extrapolation_option = ITP.ConstantExtrapolation()

# # Create the spline surrogate
itp1D = ITP.Interpolator1D(
    padding_option,
    extrapolation_option,
    buf, s, a, b; ϵ=ϵ,
) # allocates itp1D.coeffs.

# relative residual error over the interpolation region `ts`. Should be close to zero.
ts = LinRange(a, b, length(s))
q_ts = collect(ITP.query1D(u, itp1D) for u in ts)
@show norm(s - q_ts) / norm(s)

# # Visualize
# Query range.
itp_lb, itp_ub = ITP.get_itp_interval(itp1D)
border = T(2)
tq_range = LinRange(itp_lb - border, itp_ub + border, Nq)

# query.
query1D_t = collect(ITP.query1D(u, itp1D) for u in tq_range)

f_tq = f.(tq_range) # oracle

# If `Np` was set to a smaller value, then you can see some oscillations near the border.
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label="data", "x")
PLT.plot(tq_range, f_tq, label="oracle")
PLT.plot(tq_range, query1D_t, label="query1D()")
PLT.legend()
PLT.title("Spline model")

# # Quadratic Extrapolation
# One can approach non-constant extrapolation in the 1D setting as joining the spline model for querying the interpolation region with another function that handles the extrapolation region. One can use a common 1D function joining technique that ensures a differentiable transition; see [this video]( https://www.youtube.com/watch?v=vD5g8aVscUI&t=523s).

# For the 1D, real-valued setting, this package provides a quadratic extrapolation model based on the function joining technique. This is implemented by the `QuadraticExtrapolator1D` constructor, which wraps an existing spline model, in this case `itp1D`, with a quadratic model for each boundary. Let's make such a model with a transition interval equivalent to 1 sampling period between samples.
num_transition_samples = 1
etp1D = ITP.QuadraticExtrapolator1D(itp1D; num_transition_samples)
q_tq = collect(ITP.query1D(u, etp1D) for u in tq_range)

# There are two boundaries in the 1D case: the upper and lower bounds. The quadratic model of a boundary is based on the value and first and second derivatives of the border time stamps of the interpolation samples, i.e., `ts[begin]` and `ts[end]` in this case. You can also get those time stamps by doing 
@show ITP.get_itp_interval(itp1D)

# You can see the extrapolation region is the quadratic model.
PLT.figure(fig_num)
fig_num += 1
PLT.plot(t_range, s, label="data", "x")
PLT.plot(tq_range, f_tq, label="oracle")
PLT.plot(tq_range, q_tq, label="query1D()")
PLT.legend()
PLT.title("Quadratic extrapolation model")
# We can lengthen the transition interval by setting `num_transition_samples` to a larger number. Try setting it to 10, and re-run the plot.

# ## Performance
# 
# For inputs in the interpolation interval, i.e. within `get_itp_interval(itp1D)`, only the spline model needs to be evaluated when querying the quadratic extrapolation model. However, it is slower than querying the spline model due to the additional overhead for destructuring the bundled data structures and additional bounds checking. Here are the timings on the test machine:
itp_lb, itp_ub = ITP.get_itp_interval(itp1D)
x_test = (itp_lb + itp_ub) / 2
@btime ITP.query1D($x_test, $itp1D) # direct call to spline model.
@btime ITP.query1D($x_test, $etp1D) # quadratic extrapolation.

# The join transition intervals are the extrapolation intervals next to the interval interval, with an interval length specified by the `num_transition_samples` keyword from the `QuadraticExtrapolator1D` constructor.
# For inputs in a transition interval, the quadratic extrapolation model needs to query both the spline model and the qaudratic model. This makes it slower than querying the spline model, because the spline model only implements constant extrapolation without function joining.
sampling_period = ts[2] - ts[1]
x_test = itp_ub + sampling_period * num_transition_samples / 2
@btime ITP.query1D($x_test, $itp1D) # direct call to spline model.
@btime ITP.query1D($x_test, $etp1D) # quadratic extrapolation.
# For inputs that are not in a transition interval and are not in the interpolation interval, the quadratic extrapolation model only needs to query the quadratic model. Since the quadratic model does not require type casting and truncation from floating-point to integers, on the test machine, the quadratic extrapolation model is faster than the spline model.
x_test = itp_ub + num_transition_samples + 1
@btime ITP.query1D($x_test, $itp1D) # direct call to spline model.
@btime ITP.query1D($x_test, $etp1D) # quadratic extrapolation.

nothing