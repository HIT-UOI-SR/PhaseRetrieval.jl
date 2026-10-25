# PhaseRetrieval

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HIT-UOI-SR.github.io/PhaseRetrieval.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HIT-UOI-SR.github.io/PhaseRetrieval.jl/dev)
[![Build Status](https://github.com/HIT-UOI-SR/PhaseRetrieval.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/HIT-UOI-SR/PhaseRetrieval.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/HIT-UOI-SR/PhaseRetrieval.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HIT-UOI-SR/PhaseRetrieval.jl)

Classic phase retrieval algorithms.

Requires Julia 1.10 or later. The package currently provides the Gerchberg–Saxton
algorithm, with FFT propagation by default and a configurable `projector` pair.
OpticalPropagation 0.4 is supported through this interface; see the
[optical propagation example](https://HIT-UOI-SR.github.io/PhaseRetrieval.jl/dev/#Optical-propagation).

## Development

Instantiate the package environment and run the tests, including integration
tests with OpticalPropagation 0.4:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Build the documentation and its executable examples with the local package:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

CI tests Julia 1.10 and the latest stable release on Linux, macOS, and Windows,
plus an experimental nightly run on Linux. Documentation deploys on pushes to
`master` and version tags. Package, test, and documentation environments declare
their own compatibility bounds; generated manifests stay local so each supported
Julia version can resolve compatible dependencies. Dependabot checks Julia
dependencies and GitHub Actions weekly.
