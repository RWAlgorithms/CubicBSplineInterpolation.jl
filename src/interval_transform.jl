# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

struct IntervalConversion{T <: AbstractFloat}
    a::T
    d_div_bma::T
end

function IntervalConversion(a::T, b::T, N_samples::Int) where {T <: AbstractFloat}
    d = N_samples - 1
    return IntervalConversion(a, d / (b - a))
end

# not used.
"""
    to_std_interval(x::T, a::T, b::T, d::Int) where T <: AbstractFloat

converts x ∈ [a,b] to u ∈ [0,d].
If `N` is the number of samples, then `d` should be
"""
function to_std_interval(x::T, a::T, b::T, d::Int) where {T <: AbstractFloat}

    #return (x-a)*(d-c)/(b-a)+c if we're going to interval [c,d]
    return (x - a) * d / (b - a)
end

# use.
function to_std_interval(x::T, a::T, d_div_bma::T) where {T <: AbstractFloat}
    return (x - a) * d_div_bma
end

# not used
"""
    from_std_interval(x::T, a::T, b::T, d::Int) where T <: AbstractFloat

converts x ∈ [a,b] to u ∈ [0,d].
If `N` is the number of samples, then `d` should be
"""
function from_std_interval(u::T, a::T, b::T, d::Int) where {T <: AbstractFloat}

    # (x-a)*d/(b-a) = u
    return u * (b - a) / d + a
end

# use.
function from_std_interval(u::T, a::T, d_div_bma::T) where {T <: AbstractFloat}
    # u = (x-a)*d_div_bma
    return u / d_div_bma + a
end
