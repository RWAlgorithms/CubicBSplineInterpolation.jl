# CubicBSplineInterpolation.jl
Minimalist separable cubic b-spline interpolation for samples on a grid. Example scripts are in `examples/`. More documentation coming soon.

A border of 4 pixels exist where this package does not attempt to fit to the data, since it is the border region. See the plot in the `1D_demo.jl` for an illustration. This is an intensional design to avoid complicated options and branching code for handling the border regions, which slows down computational for queries. See [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for an interpolation package with more functionality.

The main methods are `Interpolator1D` to create the interpolator model, `query_interior` to query model, and `get_query_interval` to get the region where the model is fitted against the provided samples.

No `@inbounds` nor `@fastmath` is used in this code, and the focus is on minmal dependencies with only one type signatures in the composite types used by this model, and to ensure type stability. The `query_interior` is currently a little slower than the callable interpolation objects from [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl); see the notes at the end of the `example/demo_1D.jl` and `example/demo_2D.jl` for the timing on a Ryzen 7 1700 cpu in Fedora Workstation Release 40, Julia v1.11-rc1.
