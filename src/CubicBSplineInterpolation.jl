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


include("interval_transform.jl")
include("cubic1D.jl")
include("cubic2D.jl")

export IntervalConversion, get_coeffs, query1D, query2D, Interpolator1D, Interpolator2D, Interpolator1DComplex, Interpolator2DComplex

end # module CubicBSplineInterpolation
