# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module CubicBSplineInterpolation

function twothirds(::Type{T}) where T <: AbstractFloat
    return T(2)/T(3)
end

# function twothirds(::Type{Float64})
#     return 0.6666666666666666
# end

# function twothirds(::Type{Float32})
#     return 0.6666667f0
# end

function half(::Type{T}) where T <: AbstractFloat
    return one(T)/T(2)
end


include("padding.jl")

include("interval_transform.jl")

include("cubic1D.jl")
include("cubic2D.jl")
include("complex.jl") # complex-valued version.

include("derivatives1D.jl")
include("derivatives2D.jl")

export IntervalConversion, 
get_coeffs, query1D, query2D,
FitBuffer1D, SetupBuffer2D,
Interpolator1D, Interpolator2D,
Interpolator1DComplex, Interpolator2DComplex,

# padding options
LinearPadding, ConstantPadding, Lagrange4Padding,

# extrapolation options
ZeroExtrapolation, ConstantExtrapolation

# update coeffs directly
public update_coeffs!, get_num_coeffs,
query1D_parameter_derivatives

end # module CubicBSplineInterpolation
