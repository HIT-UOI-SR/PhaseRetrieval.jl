# PhaseRetrieval

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HIT-UOI-SR.github.io/PhaseRetrieval.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HIT-UOI-SR.github.io/PhaseRetrieval.jl/dev)
[![Build Status](https://github.com/HIT-UOI-SR/PhaseRetrieval.jl/workflows/CI/badge.svg)](https://github.com/HIT-UOI-SR/PhaseRetrieval.jl/actions)
[![Coverage](https://codecov.io/gh/HIT-UOI-SR/PhaseRetrieval.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HIT-UOI-SR/PhaseRetrieval.jl)

Classic phase retrieval algorithms.

| Algorithm | Function | Measurements / constraints |
| --- | --- | --- |
| Gerchberg–Saxton | `gerchbergsaxton(Ain, Aout)` | Two amplitudes, forward/inverse propagation |
| Yang–Gu (杨–顾) | `yanggu(Ain, Aout, G)` | Two amplitudes, general linear transform (including nonunitary) |
| Hybrid input–output | `hio(Aout, support)` | Diffraction amplitude, Boolean support, optional nonnegativity |
| Relaxed averaged alternating reflections | `raar(Aout, support)` | Same constraints as HIO |
| Single-beam multiple-intensity reconstruction | `sbmir(intensities, positions; wavelength, pixelsize)` | Intensities at multiple propagation planes |

GS, Yang–Gu, HIO and RAAR take **amplitudes**, not intensities. Convert camera
intensity data with `sqrt.(I)`. SBMIR takes **intensities** directly. Arrays must
contain finite, nonnegative measurements; all solvers leave their inputs intact.

```julia
using PhaseRetrieval, FFTW, Random
Random.seed!(42)

# Two-plane Fourier measurements, with an independent initial phase.
truth = (0.5 .+ rand(32, 32)) .* exp.(im .* randn(32, 32))
Ain, Aout = abs.(truth), abs.(fft(truth))
field, predicted = gerchbergsaxton(Ain, Aout; maxiter=500)
phase = angle.(field)

# Support-based diffraction imaging of a real nonnegative object.
support = falses(32, 32)
support[12:20, 12:20] .= true
object = rand(32, 32) .* support
A = abs.(fft(object))
xhio, _ = hio(A, support; nonnegative=true, beta=0.95, maxiter=500)
xraar, _ = raar(A, support; nonnegative=true, beta=0.95, maxiter=500)
estimate = ifelse.(support, max.(real.(xraar), 0), 0.0)

# General (possibly rectangular/nonunitary) matrix on vec(field).
G = randn(ComplexF64, 12, 4)
u = exp.(im .* randn(4))
uyg, vyg = yanggu(abs.(u), abs.(G*u), G; maxiter=500, inneriter=10)

# Multiple planes: units here are metres; spacings may be unequal.
z = [0.0, 0.8e-3, 2.1e-3, 3.9e-3, 6.2e-3]
optics = (wavelength=532e-9, pixelsize=(8e-6, 8e-6))
intensities = [abs2.(fresnelpropagate(truth, d; optics...)) for d in z]
fields = sbmir(intensities, z; optics..., init=ones(size(truth)), maxiter=200)
```

The first four solvers return `(field, forward(field))`. The second array is a
physical prediction, so it need not have the measured amplitude before
convergence. HIO/RAAR return the raw iterate, including feedback outside the
support. SBMIR returns a vector of fields propagated from a common first-plane
estimate. See the [algorithm and API documentation](docs/src/index.md) for
equations, custom propagators, stopping rules and references.

**Compatibility note:** GS retains its function name, keywords and tuple return,
but now copies `init`, defines `maxiter=0`, always enforces `Ain`, and returns the
actual propagation of its final field instead of the last amplitude-clamped
intermediate. Its stopping error is measured against `Aout` in the output plane.

These are nonconvex iterative methods. A small measurement residual does not
establish a unique phase, and convergence depends on data, initialization and
parameters (including RAAR's `beta`). Global phase is unobservable; support-only
Fourier data can also admit translations or conjugate inversion. Use multiple
initializations and inspect residuals for experimental data.

Run the formula, reconstruction and validation tests with:

```julia
using Pkg
Pkg.test()
```
