# CubicBSplineInterpolation.jl
Minimalist separable cubic b-spline interpolation for samples on a grid. Example scripts are in `examples/`. No derivatives are provided in the current implementation.

# License
This project is licensed under the Mozilla Public License v2.0; see the LICENSE file for details. Individual source files may contain the following tag instead of the full license text:
```
SPDX-License-Identifier: MPL-2.0
```

Using SPDX enables machine processing of license information based on the SPDX License Identifiers and makes it easier for developers to see at a glance which license they are dealing with.

# Notable attributes
The forward-backward digital filtering algorithm from Ref. [1] (see the last section of this README) was used. The derivation there is based on discrete Fourier transform, which assumes the underlying continuous function is periodic. When the provided interpolation samples are from a non-periodic continuous function, this boundary condition often yield an interpolation result that has some oscillatory behavior.

To reduce implementation complexity arising from border coefficients, the query function for *CubicBSplineInterpolation.jl* clamps the query input such that out-of-bounds does not occur on the coefficient matrix. The support of the symmetric spline implemented takes 4 sample lengths, so the query algorithm needs a padded border of at least 2 samples. The upside for this design choice is slightly faster querying times due to less logic branching.

Both the oscillatory and border issue are alleviated by having the algorithm automatically pad the user-provided interpolation samples at the boundaries. See [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for an interpolation package that is based on a more conventional spline nomenclature and more boundary conditions.

# Usage notes
- First, allocate a buffer of type `FitBuffer1D` or `FitBuffer2D`, then use it to create the interpolation coefficients, which is of type `Interpolator1D`, `Interpolator1DComplex`, `Interpolator2D`, or `Interpolator2DComplex`.

- Query the interpolation via `query1D` or `query2D`, which requires the interpolation coefficients.

- Update the coefficients via `update_itp!`.

## Sample padding length
When one creates a `FitBuffer1D` or `FitBuffer2D` variable, one can pass to their constructors the optional keyword `N_padding` to specify the padding length.

The larger `N_padding` is, the less oscillatory behavior the interpolation result should be at the boundary. It needs to be larger than or equal to `5`, because this library's query code needs two padded samples to avoid out-of-bounds access, and another three samples to transition to a extrapolation constant.

The default is `N_padding = 10` for these constructors.

## Sample padding options
`LinearPadding()`, `ConstantPadding`, `Lagrange4Padding` specifies the linear, constant, and 4-th order Lagrange extrapolation on the boundary samples to pad the user-provided samples. In order of fast to slow computational time: `ConstantPadding()`, `LinearPadding()`, `Lagrange4Padding`.

The default is `LinearPadding()`.

## Extrapolation options
`ZeroExtrapolation()` and `ConstantExtrapolation()` transitions the query function from the last fit sample to either the zero function or a constant function, respectively. The transition is part of the spline surrogate, so it is differentiable. The constant for the `ConstantExtrapolation()` option depends on the value of `N_padding` and the padding option used. 

The default is `ConstantExtrapolation()`

# Examples

## Example: 1D, Real-valued
```julia
import CubicBSplineInterpolation as ITP
using LinearAlgebra

T = Float64

ϵ = T(1e-8) # numerical tolerance for the digital filtering algorithm.

# interpolation interval
N = 1000
a = T(-10)
b = T(3.4)
t_range = LinRange(a, b, N)

# specify the oracle function.
f = xx->exp(-(1/17)*(xx/3)^2)
s = f.(t_range) # generate interpolation samples.

N_padding = 5 # This should be equal and larger than 5.
buf = ITP.FitBuffer1D(T, length(s); N_padding = N_padding)
itp1D = ITP.Interpolator1D(buf, s, a, b; ϵ = ϵ) # allocates itp1D.coeffs.

# residual error, should be near zero if N_padding >= 2.
s_rec = collect( ITP.query1D(u, itp1D) for u in t_range )
@show norm(s-s_rec)/norm(s)

# query.
Nq = 3333
query_lb, query_ub = ITP.get_itp_interval(itp1D)
extrapolation_len = T(1.23)
tq_range = LinRange(query_lb- extrapolation_len, query_ub + extrapolation_len, Nq)
results = collect( ITP.query1D(u, itp1D) for u in tq_range )
# the extrapolated result should be clamped to the border interpolation samples.

# sanity check.
ITP.update_itp!(itp1D, buf, s; ϵ = ϵ) # fits and overwrites existing coefficients.
results2 = collect( ITP.query1D(u, itp1D) for u in tq_range )
@show norm(results - results2) # should be zero.

# Allocation benchmark for allocations.
using BenchmarkTools
@btime ITP.Interpolator1D($buf, $s, $a, $b; ϵ = $ϵ) # allocates. Used only to generate `itp1D`.
@btime ITP.update_itp!($itp1D, $buf, $s; ϵ = $ϵ) # Does not allocate.

```

## Example: 1D, complex-valued
```julia
import CubicBSplineInterpolation as ITP
using LinearAlgebra

T = Float64

ϵ = T(1e-8) # numerical tolerance for the digital filtering algorithm.

# generate data.
f_real = (xx)->sinc((xx)^2)
f_imag = (xx)->tanh((xx)^2)
f = (xx)->Complex(sinc((xx)^2), tanh((xx)^2))
t_range = LinRange(T(-3), T(3), 1000)
Sr = [f_real(x1) for x1 in t_range]
Si = [f_imag(x1) for x1 in t_range]

# interpolate.
cNp = 5
cbuf = ITP.FitBuffer1D(T, length(Sr); N_padding = cNp)
citp = ITP.Interpolator1DComplex(cbuf, Sr, Si, first(t_range), last(t_range); ϵ = ϵ)

# fit residual error. Should be near zero.
qc = xx->ITP.query1D(xx, citp)
qc_s = qc.(t_range)
println("Interpolation fit residual:")
@show norm(Sr .+ im .* Si - qc_s)/norm(Sr .+ im .* Si) 

# query
query_lb, query_ub = ITP.get_itp_interval(citp)
tq_range = LinRange(query_lb, query_ub, 10000)
results = [ ITP.query1D(x1, citp) for x1 in tq_range ]


# sanity check.
ITP.update_itp!(citp, cbuf, Sr, Si; ϵ = ϵ) # fits and overwrites existing coefficients.
results2 = collect( ITP.query1D(u, citp) for u in tq_range )
@show norm(results - results2) # should be zero.

# Allocation benchmark for allocations.
using BenchmarkTools
a, b = first(t_range), last(t_range)
@btime ITP.Interpolator1DComplex($cbuf, $Sr, $Si, $a, $b; ϵ = $ϵ) # allocates. Used only to generate `citp`.
@btime ITP.update_itp!($citp, $cbuf, $Sr, $Si; ϵ = $ϵ) # Does not allocate.
```

## Example: 2D, real-valued
```julia
import CubicBSplineInterpolation as ITP
using LinearAlgebra

T = Float32
ϵ = eps(T)*2 # numerical tolerance for the digital filtering algorithm.

# generate data
a1 = T(-4.0)
b1 = T(3.45)
N1 = 100
a2 = T(-2.0)
b2 = T(1.23)
N2 = 113
t_range1 = LinRange(a1, b1, N1)
t_range2 = LinRange(a2, b2, N2)
f = (xx,yy)->sinc((xx)^2+(yy)^2)
S = [f(x1,x2) for x1 in t_range1, x2 in t_range2]

# fit coefficients.
buf = ITP.FitBuffer2D(T, size(S); N_padding = (10,9))
itp2D = ITP.Interpolator2D(buf, S, a1, b1, a2, b2; ϵ = ϵ)

# fit residual.
q_S = [ITP.query2D(x1, x2, itp2D) for x1 in t_range1, x2 in t_range2]
@show norm(S - q_S)/norm(S)

# query
Nq1 = N1*33
Nq2 = N2*33
extrapolation_len = T(0.987)
aq1, bq1 = t_range1[begin] - extrapolation_len, t_range1[end] + extrapolation_len
aq2, bq2 = t_range2[begin] - extrapolation_len, t_range2[end] + extrapolation_len
tq_range1 = LinRange(aq1, bq1, Nq1)
tq_range2 = LinRange(aq2, bq2, Nq2)

results = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]

# sanity check.
ITP.update_itp!(itp2D, buf, S; ϵ = ϵ) # fits and overwrites existing coefficients.
results2 = [ITP.query2D(x1, x2, itp2D) for x1 in tq_range1, x2 in tq_range2]
@show norm(results - results2) # should be zero.

# Allocation benchmark for allocations.
using BenchmarkTools
@btime ITP.Interpolator2D($buf, $S, $a1, $b1, $a2, $b2; ϵ = $ϵ) # allocates. Used only to generate `itp2D`.
@btime ITP.update_itp!($itp2D, $buf, $S; ϵ = $ϵ) # Does not allocate.
```

## Example: 2D, complex-valued
```julia
import CubicBSplineInterpolation as ITP
using LinearAlgebra

T = Float32
ϵ = eps(T)*2 # numerical tolerance for the digital filtering algorithm.

# specify oracle.
f_real = (xx,yy)->sinc((xx)^2+(yy)^2)
f_imag = (xx,yy)->tanh((xx)^2+(yy)^2)
f = (xx,yy)->Complex(sinc((xx)^2+(yy)^2), tanh((xx)^2+(yy)^2)) # This allocates for some reason. Perhaps due to "boxing". Complex(f_real(xx,yy), f_imag(xx,yy))

t_range1 = LinRange(T(-3), T(3), 1000)
t_range2 = LinRange(T(-1.23), T(4.56), 783)

# Generate samples.
Sr = [f_real(x1,x2) for x1 in t_range1, x2 in t_range2]
Si = [f_imag(x1,x2) for x1 in t_range1, x2 in t_range2]

# fit coefficients.
cbuf = ITP.FitBuffer2D(T, size(Sr); N_padding = (8,7))
citp = ITP.Interpolator2DComplex(ITP.LinearPadding(), ITP.ConstantExtrapolation(), cbuf, Sr, Si, first(t_range1), last(t_range1), first(t_range2), last(t_range2))

# This also works: ITP.Interpolator2DComplex(cbuf, Sr, Si, first(t_range1), last(t_range1), first(t_range2), last(t_range2))
# Without specifying the type of padding strategy, it defaults to LinearPadding.

# fit residual error. Should be near zero.
qc_s = [ITP.query2D(x1, x2, citp) for x1 in t_range1, x2 in t_range2]
println("Interpolation fit residual:")
@show norm(Sr .+ im .* Si - qc_s)/norm(Sr .+ im .* Si) 

# query.
query_lb1, query_ub1, query_lb2, query_ub2 = ITP.get_itp_interval(citp)
extrapolation_len = T(0.456)
tq_range1 = LinRange(query_lb1 - extrapolation_len, query_ub1 + extrapolation_len, 10000)
tq_range2 = LinRange(query_lb2 - extrapolation_len, query_ub2 + extrapolation_len, 1473)

results = [ ITP.query2D(x1, x2, citp) for x1 in tq_range1, x2 in tq_range2 ]

# sanity check.
ITP.update_itp!(citp, cbuf, Sr, Si; ϵ = ϵ) # fits and overwrites existing coefficients.
results2 = [ITP.query2D(x1, x2, citp) for x1 in tq_range1, x2 in tq_range2]
@show norm(results - results2) # should be zero.

# Allocation benchmark for allocations.
using BenchmarkTools
a1, a2, b1, b2 = first(t_range1), last(t_range1), first(t_range2), last(t_range2)
@btime ITP.Interpolator2DComplex($cbuf, $Sr, $Si, $a1, $b1, $a2, $b2; ϵ = $ϵ) # allocates. Used only to generate `citp`.
@btime ITP.update_itp!($citp, $cbuf, $Sr, $Si; ϵ = $ϵ) # Does not allocate.
```

# Derivatives
The cubic B-spline has continuous first and second derivatives, but discontinuous third derivatives.

# ## First derivatives
If `itp` is a variable of type or `Interpolator1D`, then use the following to get the derivative, `g`:
```julia
g = query1D_derivative1(x, itp)
```

If `itp` is a variable of type or `Interpolator1DComplex`, then use the following to get the real part's derivative `gr` and the imaginary part's derivaitve `gi`:
```julia
gr, gi = query1D_derivative1(x, itp)
```

If `itp` is a variable of type or `Interpolator2D`, then use the following to get the derivative in the first and second coordinate axis directions, stored in `g1`, `g2`, respectively:
```julia
g1, g2 = query2D_derivative1(x1, x2, itp)
```

If `itp` is a variable of type or `Interpolator2DComplex`, then use the following to get the real part's derivative `gr` and the imaginary part's derivaitve `gi`:
```julia
g1r, g2r, g1i, g2i = query2D_derivative1(x1, x2, itp)
```
The last letter for the left-hand side variables indicate whether the variable is for the real (`r`) or imaginary (`i`) part.

# ## Second derivatives
If `itp` is a variable of type or `Interpolator1D`, then use the following to get the second derivative, `g2`:
```julia
g2 = query1D_derivative2(x, itp)
```

If `itp` is a variable of type or `Interpolator1DComplex`, then use the following to get the real part's second derivative `gr` and the imaginary part's second derivaitve `gi2`:
```julia
gr2, gi2 = query1D_derivative2(x, itp)
```

If `itp` is a variable of type or `Interpolator2D`, then use the following to get the entries of the Hessian matrix:
```julia
g11, g12, g21 = query2D_derivative2(x1, x2, itp)
```
where the Hessian matrix is:
```
H = [
    g11 g12;
    g21 g22;
]
```
Schwarz's (or Clairaut's) theorem on equality of mixed partials mean `g12 = g21`.

If `itp` is a variable of type or `Interpolator2DComplex`, then use the following to get the real part's derivative `gr` and the imaginary part's derivaitve `gi`:
```julia
g11r, g12r, g22r, g11i, g12i, g22i = query2D_derivative2(x1, x2, itp)
```
The last letter for the left-hand side variables indicate whether the variable is for the real (`r`) or imaginary (`i`) part.

# Work in progress
The example scripts in `/examples` are being revised. So far, `demo1D` has been processed. The original 1D and 2D example scripts can be found in `/examples/AD/`.

# Reference
1. Michael Unser, Splines: a perfect fit for signal and image processing, 1999, DOI: 10.1109/79.799930.
