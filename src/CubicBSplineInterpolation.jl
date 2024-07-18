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
include("cubic.jl")

export IntervalConversion, get_coeffs, query_interior, Interpolator1D, Interpolator2D

end # module CubicBSplineInterpolation
