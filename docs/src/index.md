```@meta
CurrentModule = PhaseRetrieval
```

# PhaseRetrieval

Classic phase retrieval algorithms.

## Optical propagation

On Julia 1.10 or later, install PhaseRetrieval, OpticalPropagation 0.4, and Unitful
in your environment. Pass a forward/backward propagation pair to the existing
`projector` keyword. OpticalPropagation's `fresnel2` uses the same physical grid
in both directions, and its curried form accepts unitful propagation distances.

The following executable example starts with a zero phase estimate. The input
is a `MonoLightField2D`, so wavelength, window size, and window centre are retained
during the iterations:

```@example optical
using PhaseRetrieval, OpticalPropagation, Unitful

data = [exp(-3 * (x^2 + y^2)) * cis(0.8x^2 - 0.6y)
        for x in range(-1, 1; length=9), y in range(-1, 1; length=8)]
truth = MonoLightField2D(data; wavelength=532u"nm",
                        size=(0.8u"mm", 0.6u"mm"),
                        center=(0.2u"mm", -0.1u"mm"))
forward, backward = fresnel2(20u"mm"), fresnel2(-20u"mm")
measured = forward(truth)
recovered, output = gerchbergsaxton(abs.(truth), abs.(measured);
    projector=(forward, backward),
    init=MonoLightField2D(truth; data=abs.(parent(truth))),
    maxiter=30)

@assert physicalsize(recovered) == physicalsize(truth)
@assert gridcenter(recovered) == gridcenter(truth)
@assert abs.(parent(recovered)) ≈ abs.(parent(truth))
nothing # hide
```

Use matching sampling grids for the amplitude measurements and propagators.
The two returned fields satisfy their respective amplitude constraints; at a
finite iteration count they need not be an exact forward-propagated pair.
OpticalPropagation and Unitful are test/example dependencies, so the default FFT
solver only loads FFTW and Distances.

## API

```@index
```

```@autodocs
Modules = [PhaseRetrieval]
```
