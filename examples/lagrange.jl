# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# test out non-allocating Lagrange extrapolation.

using Random

PLT.close("all")
fig_num = 1

T = Float64

s = Memory{T}([0.9613874912649977, 0.9729459646121265, 0.9825234242313965, 0.9900589251576111, 0.9955042335906311, 0.9988243410766205, 0.9999978393584623, 0.9990171502802508, 0.9958886072917987, 0.9906323873268043])
ts = LinRange(T(2.1818181818181817), T(3.4), 10)

xs = ts[end-4:end]
ys = s[end-4:end]
# randn!(Random.Xoshiro(0), ys)



tqs = LinRange(T(2), T(4), 100)

q = xx->ITP.lagrange4_etp(xx, xs, ys)
q_tqs = q.(tqs)

PLT.figure(fig_num)
fig_num += 1
PLT.plot(xs, ys, "o", label = "data")
PLT.plot(tqs, q_tqs, label = "etp")
PLT.title("3rd-order Lagrange extrapolation.")

# # 

e1, e2, e3, e4 = ITP.eval_reflection_etp(ts, s)


PLT.figure(fig_num)
fig_num += 1
PLT.plot([e1; e2; s; e3; e4], "o")
PLT.title("extended data, append 2 on either boundary")

nothing