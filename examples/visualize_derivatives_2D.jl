# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

if !isdefined(Main, :CubicBSplineInterpolation)
    include("a.jl")
end

PLT.close("all")
fig_num = 1

#const T = Float64
T = Float32

ϵ = eps(T)*2
a1 = T(-4.0)
b1 = T(3.45)
N1 = 100
a2 = T(-2.0)
b2 = T(1.23)
N2 = 113
t_range1 = LinRange(a1, b1, N1)
t_range2 = LinRange(a2, b2, N2)
f = (xx,yy)->sinc((xx)^2+(yy)^2+20)
S = [f(x1,x2) for x1 in t_range1, x2 in t_range2]

# fit.
buf = ITP.FitBuffer2D(T, size(S); N_padding = (10,10))
#itp2D = ITP.Interpolator2D(buf, S, a1, b1, a2, b2; ϵ = ϵ)

padding_option = ITP.LinearPadding()
extrapolation_option = ITP.ZeroExtrapolation()

padding_option = ITP.LinearPadding()
extrapolation_option = ITP.ConstantExtrapolation()

itp2D = ITP.Interpolator2D(
    padding_option,
    extrapolation_option,
    buf, S, a1, b1, a2, b2; ϵ = ϵ,
)

# query.
Nq1 = N1*10
Nq2 = N2*10
aq1, bq1 = t_range1[begin], t_range1[end]
aq2, bq2 = t_range2[begin], t_range2[end]

extrapolation_len = T(0.3)
tq_range1 = LinRange(aq1 - extrapolation_len, bq1 + extrapolation_len, Nq1)
tq_range2 = LinRange(aq2 - extrapolation_len, bq2 + extrapolation_len, Nq2)

Yq = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq)
PLT.title("extrapolation length: $(extrapolation_len), Yq")

# ## First derivative

dq_tq = [ ITP.query2D_derivative1(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2 ]

h = xx->ITP.query2D(xx[1], xx[2], itp2D)
dq_tq_ND = [ FiniteDiff.finite_difference_gradient(h, [x1; x2]) for x1 in tq_range1, x2 in tq_range2 ]

dq_tq_x1 = map(xx->xx[begin], dq_tq)
dq_tq_ND_x1 = map(xx->xx[begin], dq_tq_ND)
dq_tq_x2 = map(xx->xx[begin+1], dq_tq)
dq_tq_ND_x2 = map(xx->xx[begin+1], dq_tq_ND)


println("First derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(dq_tq_x1 - dq_tq_ND_x1)/norm(dq_tq_ND_x1)
@show norm(dq_tq_x2 - dq_tq_ND_x2)/norm(dq_tq_ND_x2)

dq_tq_vec = [ [u[1]; u[2]] for u in dq_tq ]
PLT.figure(fig_num)
fig_num += 1
PLT.imshow(norm.(dq_tq_ND .- dq_tq_vec))
PLT.colorbar()
PLT.title("norm difference between analytical and numerical derivatives: gradient")

# ### Visualize slice at the largest norm discrepancy, to see if the discrepancy is plausible.
_, coordinate = findmax(norm.(dq_tq_ND .- dq_tq_vec))

x1_fixed, x2_fixed = tq_range1[coordinate[1]], tq_range2[coordinate[2]]
mq1 = LinRange(first(tq_range1), last(tq_range1), 10_000)
mq2 = LinRange(first(tq_range2), last(tq_range2), 10_000)

dq_tq_1 = [ ITP.query2D_derivative1(x1, x2_fixed, itp2D)[1] for x1 in mq1 ]
dq_tq_ND_1 = [ FiniteDiff.finite_difference_gradient(h, [x1; x2_fixed])[1] for x1 in mq1 ]

dq_tq_2 = [ ITP.query2D_derivative1(x1_fixed, x2, itp2D)[end] for x2 in mq2 ]
dq_tq_ND_2 = [ FiniteDiff.finite_difference_gradient(h, [x1_fixed; x2])[2] for x2 in mq2 ]


h1 = xx->h([xx; x2_fixed])
h2 = xx->h([x1_fixed; xx])

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, h1.(mq1))
PLT.xlabel("x1")
PLT.title("the query function, slice at x2 = $(x2_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, dq_tq_1, label = "implemented")
PLT.plot(mq1, dq_tq_1, "x")
PLT.plot(mq1, dq_tq_ND_1, "--", label = "numerical")
PLT.xlabel("x1")
PLT.title("1st derivative: numerical vs analytical, slice at x2 = $(x2_fixed)")
PLT.legend()

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, h2.(mq2))
PLT.xlabel("x2")
PLT.title("the query function, slice at x1 = $(x1_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, dq_tq_2, label = "implemented")
PLT.plot(mq2, dq_tq_2, "x")
PLT.plot(mq2, dq_tq_ND_2, "--", label = "numerical")
PLT.xlabel("x2")
PLT.title("1st derivative: numerical vs analytical, slice at x1 = $(x1_fixed)")
PLT.legend()

"""
The implemetned analytical gradient is continuous and roughly matches the numerical, so it is plausible the implementation is correct.
"""


# ## Second derivative

d2q_tq = [ ITP.query2D_derivative2(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2 ]

h = xx->ITP.query2D(xx[1], xx[2], itp2D)
d2q_tq_ND = [ FiniteDiff.finite_difference_hessian(h, [x1; x2]) for x1 in tq_range1, x2 in tq_range2 ]

d2q_tq_x1x1 = map(xx->xx[begin], d2q_tq)
d2q_tq_ND_x1x1 = map(xx->xx[1,1], d2q_tq_ND)

d2q_tq_x1x2 = map(xx->xx[begin+1], d2q_tq)
d2q_tq_ND_x1x2 = map(xx->xx[1,2], d2q_tq_ND)

d2q_tq_x2x2 = map(xx->xx[begin+2], d2q_tq)
d2q_tq_ND_x2x2 = map(xx->xx[2,2], d2q_tq_ND)

println("First derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(d2q_tq_x1x1 - d2q_tq_ND_x1x1)/norm(d2q_tq_ND_x1x1)
@show norm(d2q_tq_x1x2 - d2q_tq_ND_x1x2)/norm(d2q_tq_ND_x1x2)
@show norm(d2q_tq_x2x2 - d2q_tq_ND_x2x2)/norm(d2q_tq_ND_x2x2)

#@assert 1==23

d2q_tq_vec = [ [u[1] u[2]; u[2] u[3]] for u in d2q_tq ]
PLT.figure(fig_num)
fig_num += 1
PLT.imshow(norm.(d2q_tq_ND .- d2q_tq_vec))
PLT.colorbar()
PLT.title("norm difference between analytical and numerical derivatives: hessian")


# ## Visualize slice at the largest norm discrepancy, to see if the discrepancy is plausible.
# ### slice along x2 direction.
_, coordinate = findmax(norm.(d2q_tq_ND .- d2q_tq_vec))

x1_fixed, x2_fixed = tq_range1[coordinate[1]], tq_range2[coordinate[2]]
mq1 = LinRange(first(tq_range1), last(tq_range1), 10_000)
mq2 = LinRange(first(tq_range2), last(tq_range2), 10_000)

h = xx->ITP.query2D(xx[1], xx[2], itp2D)
d2q_tq_1 = [ ITP.query2D_derivative2(x1, x2_fixed, itp2D)[1] for x1 in mq1 ]
d2q_tq_ND_1 = [ FiniteDiff.finite_difference_hessian(h, [x1; x2_fixed])[1,1] for x1 in mq1 ]

# ## need mq_cross to be a diagonal line in the input space, instead of along some axis.
# d2q_tq_12 = [ ITP.query2D_derivative2(x1, x2_fixed, itp2D)[2] for x1 in mq1 ]
# d2q_tq_ND_12 = [ FiniteDiff.finite_difference_hessian(h, [x1; x2_fixed])[1,2] for x1 in mq1 ]

d2q_tq_2 = [ ITP.query2D_derivative2(x1_fixed, x2, itp2D)[end] for x2 in mq2 ]
d2q_tq_ND_2 = [ FiniteDiff.finite_difference_hessian(h, [x1_fixed; x2])[2,2] for x2 in mq2 ]

# Schwartz symmetry theorem check.
d2q_tq_ND = [ FiniteDiff.finite_difference_hessian(h, [x1_fixed; x2]) for x2 in mq2 ]
@assert all(issymmetric.(d2q_tq_ND))

h1 = xx->h([xx; x2_fixed])
h2 = xx->h([x1_fixed; xx])

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, h1.(mq1))
PLT.xlabel("x1")
PLT.title("the query function, slice at x2 = $(x2_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, d2q_tq_1, label = "implemented")
PLT.plot(mq1, d2q_tq_1, "x")
PLT.plot(mq1, d2q_tq_ND_1, "--", label = "numerical")
PLT.xlabel("x1")
PLT.title("2nd derivative: numerical vs analytical, slice at x2 = $(x2_fixed)")
PLT.legend()

## need mq_cross to be a diagonal line in the input space, instead of along some axis.
# PLT.figure(fig_num)
# fig_num += 1
# PLT.plot(mq1, d2q_tq_1, label = "implemented")
# PLT.plot(mq1, d2q_tq_1, "x")
# PLT.plot(mq1, d2q_tq_ND_1, "--", label = "numerical")
# PLT.xlabel("x1")
# PLT.title("cross derivative: numerical vs analytical, slice at x2 = $(x2_fixed)")
# PLT.legend()

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, h2.(mq2))
PLT.xlabel("x2")
PLT.title("the query function, slice at x1 = $(x1_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, d2q_tq_2, label = "implemented")
PLT.plot(mq2, d2q_tq_2, "x")
PLT.plot(mq2, d2q_tq_ND_2, "--", label = "numerical")
PLT.xlabel("x2")
PLT.title("2nd derivative: numerical vs analytical, slice at x1 = $(x1_fixed)")
PLT.legend()

"""
The implemented 2nd derivative is continuous but itself has discontinuous derivatives.
This causes large deviations with the numerical derivatives, which are smoothed.

The 3rd derivative is discontious because Cubic B-splines are in C^{2} but not in C^{3}.
"""

# # Complex-valued.
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

# padding_option = ITP.LinearPadding()
# extrapolation_option = ITP.ZeroExtrapolation()

padding_option = ITP.LinearPadding()
extrapolation_option = ITP.ConstantExtrapolation()

cbuf = ITP.FitBuffer2D(T, size(Sr); N_padding = (8,7))
citp = ITP.Interpolator2DComplex(
    padding_option,
    extrapolation_option,
    cbuf, Sr, Si, first(t_range1), last(t_range1), first(t_range2), last(t_range2),
)

# query
extrapolation_len = T(0.5)
a1, b1, a2, b2 = ITP.get_itp_interval(citp)
tq_range1 = LinRange(a1-extrapolation_len, b1+extrapolation_len, 2000)
tq_range2 = LinRange(a2-extrapolation_len, b2+extrapolation_len, 2000)

Yq_r = [ real(ITP.query2D(x1, x2, citp)) for x1 in tq_range1, x2 in tq_range2]
Yq_i = [ real(ITP.query2D(x1, x2, citp)) for x1 in tq_range1, x2 in tq_range2]

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq_r)
PLT.colorbar()
PLT.title("real part: extrapolation length: $(extrapolation_len), Yq")


PLT.figure(fig_num)
fig_num += 1
PLT.imshow(Yq_i)
PLT.colorbar()
PLT.title("imag part: extrapolation length: $(extrapolation_len), Yq")


# ## First derivative

dq_tq = [ ITP.query2D_derivative1(x1, x2, citp) for x1 in tq_range1, x2 in tq_range2 ]
dq_tq_r = map(xx->[xx[1]; xx[2]], dq_tq)
dq_tq_i = map(xx->[xx[3]; xx[4]], dq_tq)

hr = xx->real(ITP.query2D(xx[1], xx[2], citp))
dq_tq_ND_r = [ FiniteDiff.finite_difference_gradient(hr, [x1; x2]) for x1 in tq_range1, x2 in tq_range2 ]

hi = xx->imag(ITP.query2D(xx[1], xx[2], citp))
dq_tq_ND_i = [ FiniteDiff.finite_difference_gradient(hi, [x1; x2]) for x1 in tq_range1, x2 in tq_range2 ]


println("First derivative discrepancy: Numerical and implemented analytical derivative")
@show norm(norm.(dq_tq_r .- dq_tq_ND_r))/norm(norm.(dq_tq_ND_r))
@show norm(norm.(dq_tq_i .- dq_tq_ND_i))/norm(norm.(dq_tq_ND_i))

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(norm.(dq_tq_r .- dq_tq_ND_r))
PLT.colorbar()
PLT.title("real part: derivative discrepancy between numerical and the implemented analytical")

PLT.figure(fig_num)
fig_num += 1
PLT.imshow(norm.(dq_tq_i .- dq_tq_ND_i))
PLT.colorbar()
PLT.title("imag part: derivative discrepancy between numerical and the implemented analytical")

include("slice_derivatives_2D_complex.jl")

# #### imagiary part



# TODO test the second derivative case.

nothing