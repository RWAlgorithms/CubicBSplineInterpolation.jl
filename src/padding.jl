# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# routines for generating the padded interpolation samples.

function lagrange4_etp(x::T, xs, ys::AbstractVector{T}, w0, w1, w2, w3, w4) where T <: AbstractFloat

    y0, y1, y2, y3, y4 = ys
    x0, x1, x2, x3, x4 = xs

    out_numerator = y0*w0/(x-x0) + y1*w1/(x-x1) + y2*w2/(x-x2) + y3*w3/(x-x3) + y4*w4/(x-x4)
    out_denominator = w0/(x-x0) + w1/(x-x1) + w2/(x-x2) + w3/(x-x3) + w4/(x-x4)
    return out_numerator/out_denominator
end

function lagrange4_etp(x, xs, ys)
    w0, w1, w2, w3, w4 = setup_lagrange4(xs)
    return lagrange4_etp(x, xs, ys, w0, w1, w2, w3, w4)
end

function setup_lagrange4(xs::Union{AbstractVector{T}, AbstractRange{T}}) where T <: AbstractFloat
    x0, x1, x2, x3, x4 = xs

    w0_denominator = (x0-x1)*(x0-x2)*(x0-x3)*(x0-x4)
    w0 = one(T)/w0_denominator

    w1_denominator = (x1-x0)*(x1-x2)*(x1-x3)*(x1-x4)
    w1 = one(T)/w1_denominator

    w2_denominator = (x2-x0)*(x2-x1)*(x2-x3)*(x2-x4)
    w2 = one(T)/w2_denominator

    w3_denominator = (x3-x0)*(x3-x1)*(x3-x2)*(x3-x4)
    w3 = one(T)/w3_denominator

    w4_denominator = (x4-x0)*(x4-x1)*(x4-x2)*(x4-x3)
    w4 = one(T)/w4_denominator

    return w0, w1, w2, w3, w4
end


#### boundary padding options

abstract type PaddingOption end
struct ConstantPadding <: PaddingOption end
struct Lagrange4Padding <: PaddingOption end
struct LinearPadding <: PaddingOption end

# mutates yq, output.
function _lagrange4_end_etp!(yq::AbstractVector{T}, ts::LinRange, samples::AbstractVector{T}) where T <: AbstractFloat

    Δt = step(ts)
    order = 4

    # # the larger boundary
    xs = ts[(length(ts)-order-1):end]
    ys = view(samples, (length(ts)-order-1):length(ts))
    w0, w1, w2, w3, w4 = setup_lagrange4(xs)
    
    x_offset = ts[end]
    for m in eachindex(yq)
        yq[m] = lagrange4_etp(x_offset + m*Δt, xs, ys, w0, w1, w2, w3, w4)
    end
    return nothing
end


# mutates yq, output.
# yq must be 1-indexing with stride 1.
function _lagrange4_st_etp!(yq::AbstractVector{T}, ts::LinRange, samples::AbstractVector{T}) where T <: AbstractFloat

    Δt = step(ts)
    order = 4

    # # the larger boundary
    xs = ts[begin:(order+1)]
    ys = view(samples, 1:(order+1))
    w0, w1, w2, w3, w4 = setup_lagrange4(xs)
    
    x_offset = ts[begin]
    # reverse iterator. Iterators.reverse is lazy.
    k = 0
    for m in length(yq):-1:1
        k += 1
        yq[m] = lagrange4_etp(x_offset - k*Δt, xs, ys, w0, w1, w2, w3, w4)
    end
    return nothing
end

# clamp to constant.
function _pad_samples!(::ConstantPadding, s_ext::Memory, s, Np::Integer, args...)
    st_ind = Np + 1
    fin_ind = length(s_ext) - Np

    for i in Iterators.take(eachindex(s_ext), st_ind-1)
        s_ext[i] = s[begin]
    end

    for i in (fin_ind+1):length(s_ext)
        s_ext[i] = s[end]
    end
    return nothing
end

function _pad_samples!(::Lagrange4Padding, s_ext::Memory, s, Np::Integer, xs::LinRange)
    # large boundary.
    yq = view(s_ext, (length(s_ext)-Np+1):length(s_ext))
    _lagrange4_end_etp!(yq, xs, s)

    # small boundary
    yq = view(s_ext, 1:Np)
    _lagrange4_st_etp!(yq, xs, s)
    return nothing
end

function _pad_samples!(::LinearPadding, s_ext::Memory, s, Np::Integer, args...)

    Δs = s[begin+1] - s[begin]
    k = 0
    for i in Np:-1:1
        k += 1
        s_ext[i] = s[begin] - k*Δs
    end

    Δs = s[end] - s[end-1]
    k = 0
    for i in (length(s_ext)-Np+1):length(s_ext)
        k += 1
        s_ext[i] = s[end] + k*Δs
    end
    return nothing
end


#### extrapolation
abstract type ExtrapolationOption end
struct ConstantExtrapolation <: ExtrapolationOption end
struct ZeroExtrapolation <: ExtrapolationOption end

