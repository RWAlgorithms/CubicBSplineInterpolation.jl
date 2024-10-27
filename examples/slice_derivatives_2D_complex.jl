# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

PLT.close("all")

# ### Visualize slice at the largest norm discrepancy, to see if the discrepancy is plausible.
_, coordinate = findmax(
    norm.(dq_tq_r .- dq_tq_ND_r) + norm.(dq_tq_i .- dq_tq_ND_i)
)
x1_fixed, x2_fixed = tq_range1[coordinate[1]], tq_range2[coordinate[2]]
mq1 = LinRange(first(tq_range1), last(tq_range1), 10_000)
mq2 = LinRange(first(tq_range2), last(tq_range2), 10_000)

# #### real part
dq_tq_1 = [ ITP.query2D_derivative1(x1, x2_fixed, citp)[1] for x1 in mq1 ]
dq_tq_ND_1 = [ FiniteDiff.finite_difference_gradient(hr, [x1; x2_fixed])[1] for x1 in mq1 ]

dq_tq_2 = [ ITP.query2D_derivative1(x1_fixed, x2, citp)[2] for x2 in mq2 ]
dq_tq_ND_2 = [ FiniteDiff.finite_difference_gradient(hr, [x1_fixed; x2])[2] for x2 in mq2 ]

h1 = xx->hr([xx; x2_fixed])
h2 = xx->hr([x1_fixed; xx])

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, h1.(mq1))
PLT.xlabel("x1")
PLT.title("real part: the query function, slice at x2 = $(x2_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, dq_tq_1, label = "implemented")
PLT.plot(mq1, dq_tq_1, "x")
PLT.plot(mq1, dq_tq_ND_1, "--", label = "numerical")
PLT.xlabel("x1")
PLT.title("real part: 1st derivative: numerical vs analytical, slice at x2 = $(x2_fixed)")
PLT.legend()

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, h2.(mq2))
PLT.xlabel("x2")
PLT.title("real part: the query function, slice at x1 = $(x1_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, dq_tq_2, label = "implemented")
PLT.plot(mq2, dq_tq_2, "x")
PLT.plot(mq2, dq_tq_ND_2, "--", label = "numerical")
PLT.xlabel("x2")
PLT.title("real part: 1st derivative: numerical vs analytical, slice at x1 = $(x1_fixed)")
PLT.legend()

# #### imaginary part
dq_tq_1 = [ ITP.query2D_derivative1(x1, x2_fixed, citp)[3] for x1 in mq1 ]
dq_tq_ND_1 = [ FiniteDiff.finite_difference_gradient(hi, [x1; x2_fixed])[1] for x1 in mq1 ]

dq_tq_2 = [ ITP.query2D_derivative1(x1_fixed, x2, citp)[4] for x2 in mq2 ]
dq_tq_ND_2 = [ FiniteDiff.finite_difference_gradient(hi, [x1_fixed; x2])[2] for x2 in mq2 ]

h1 = xx->hi([xx; x2_fixed])
h2 = xx->hi([x1_fixed; xx])

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, h1.(mq1))
PLT.xlabel("x1")
PLT.title("imag part: the query function, slice at x2 = $(x2_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq1, dq_tq_1, label = "implemented")
PLT.plot(mq1, dq_tq_1, "x")
PLT.plot(mq1, dq_tq_ND_1, "--", label = "numerical")
PLT.xlabel("x1")
PLT.title("imag part: 1st derivative: numerical vs analytical, slice at x2 = $(x2_fixed)")
PLT.legend()

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, h2.(mq2))
PLT.xlabel("x2")
PLT.title("imag part: the query function, slice at x1 = $(x1_fixed)")

PLT.figure(fig_num)
fig_num += 1
PLT.plot(mq2, dq_tq_2, label = "implemented")
PLT.plot(mq2, dq_tq_2, "x")
PLT.plot(mq2, dq_tq_ND_2, "--", label = "numerical")
PLT.xlabel("x2")
PLT.title("imag part: 1st derivative: numerical vs analytical, slice at x1 = $(x1_fixed)")
PLT.legend()

"""
The implemented first derivative is continuous.
When zoomed-in, looks like the second derivative is continuous.

The numerical derivative looks to be a smoothed version of the implemented analytical derivative.
"""

nothing