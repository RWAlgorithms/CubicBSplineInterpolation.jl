# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module CubicBSplineInterpolation

function twothirds(::Type{T}) where T <: AbstractFloat
    return 0.6666666666666666
end

function twothirds(::Type{Float64})
    return 0.6666666666666666
end

function twothirds(::Type{Float32})
    return 0.6666667f0
end

include("padding.jl")

include("interval_transform.jl")

include("cubic1D.jl")
include("cubic2D.jl")
include("complex.jl") # complex-valued version.


export IntervalConversion, 
get_coeffs, query1D, query2D,
FitBuffer1D, SetupBuffer2D,
Interpolator1D, Interpolator2D,
Interpolator1DComplex, Interpolator2DComplex,

# padding options
LinearPadding, ConstantPadding, Lagrange4Padding,

# extrapolation options
ZeroExtrapolation, ConstantExtrapolation

end # module CubicBSplineInterpolation
