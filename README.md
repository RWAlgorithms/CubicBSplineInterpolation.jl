# CubicBSplineInterpolation.jl
Minimalist separable cubic b-spline interpolation for samples on a grid. Example scripts are in `examples/`. More documentation coming soon.

A border of 4 pixels exist where this package does not attempt to fit to the data, since it is the border region. See the plot in the `1D_demo.jl` for an illustration. This is an intensional design to avoid complicated options and branching code for handling the border regions, which slows down computational for queries. See [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for an interpolation package with more functionality.

# Main methods
The main exported methods are:

- `Interpolator1D`, Interpolator2DComplex`, Interpolator2D`, and Interpolator2DComplex` to create the interpolator model.

- `query1D` and `query2D` to query model,

- `get_query_interval` to get the region where the model is fitted against the provided samples. Currently, it marks the region within 8 samples from the boundaries as *not* in the query interval. The choice of 8 is a heuristic, based on 4 samples being used to define the cubic spline interpolation window.

While this package offers a `get_query_interval` to give a heuristic discard of 8 samples wide border, you should do your own error assessment if low predictive error is important for your usage.

# Notable features
This package is meant for fast querying in regions away from the boundaries. For other spline interpolation use settings, please consider [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).

The interpolation of complex-valued 1D and 2D samples is done in the same function, leading to minor speed ups. These use the `Interpolator1DComplex` and `Interpolator2DComplex` composite types; see the example scripts for usage.

No `@inbounds` nor `@fastmath` is used in this code, and the focus is on minmal dependencies with only one type signatures in the composite types used by this model, and to ensure type stability.

# Error
For the limited syntehtic functions tested, it seems that the [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) has lower predictive error for query regions near the boundary. See `examples/demo_1D.jl`. This might be due to our use of a recursive filtering operation to get the interpolation coefficients. We used the procedure from Box 2 and 3 from [(Unser 1999, DOI: 10.1109/79.799930)](https://doi.org/10.1109/79.799930) without enforcing the mirror indexing on the coefficients during the convolution formula 1.

# Timing
For `Float64` samples in my tests, the `query1D` and `query2D` is currently a little slower for real-valued samples the interpolants from [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl), but a little faster for complex-valued samples. See the notes at the end of the `example/demo_1D.jl` and `example/demo_2D.jl` for the timing on the test system: Ryzen 7 1700 cpu in Fedora Workstation Release 40, Julia v1.11-rc1.

The `Float32` tests indicate `query1D` and `query2D` are a bit faster for all cases than [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) on the test system. 
