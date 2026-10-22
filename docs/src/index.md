```@meta
CurrentModule = PhaseRetrieval
```

# PhaseRetrieval

Classic phase retrieval algorithms.

## Data and return conventions

`gerchbergsaxton`, `yanggu`, `hio` and `raar` accept amplitudes (`sqrt.(I)` for
camera intensities). `sbmir` accepts intensities and takes the square roots
internally. Initialization is a field, not an array of phase angles: pass
`exp.(im .* phase)` to initialize from a phase map. Inputs are not modified.

Two-plane solvers return `(field, forward(field))`. GS and Yang–Gu enforce the
input amplitude on `field`; HIO and RAAR return their raw feedback iterate.
Their output-plane fields are predictions, not amplitude-clamped intermediates.
SBMIR returns a vector of predictions from its final first-plane field.

For an object constrained to a support, extract an estimate with
`ifelse.(support, field, zero(eltype(field)))`. If nonnegative, use
`ifelse.(support, max.(real.(field), 0), 0.0)`. Check the propagation residual of
this estimate as well as the raw iterate. A RAAR fixed point for inconsistent
constraints need not be a feasible object.

All solvers support `maxiter=0`. Setting `tol=0` forces exactly `maxiter` updates;
otherwise stopping uses the supplied `errorfunc` (mean squared difference by
default). GS and Yang–Gu check predicted output amplitude error. HIO and RAAR
also check the distance to the object constraint. SBMIR checks amplitude error
at every plane using **unconstrained forward predictions of the same field**.
Tolerances are absolute, in the units of the selected error function.

The original GS keywords and tuple structure remain supported, but the second
return value now equals the forward transform of the returned first value.
Previously it was an amplitude-clamped intermediate from a different point in
the iteration. GS now copies `init` and always enforces the input amplitude,
including early exits; its error test is now in the output plane.

## Update equations

Write ``P_a(z)=a\,\mathrm{phase}(z)``, with phase zero chosen when ``z=0``.
For a scaled unitary forward operator ``F``, let
``P_M(x)=F^{-1}P_b(Fx)``. The default `(fft, ifft)` uses Julia's unshifted FFT
ordering, with no implicit `fftshift` or intensity normalization.

**Gerchberg–Saxton:** ``x_{k+1}=P_a(P_M(x_k))``. A custom non-mutating pair is
passed as `projector=(forward, inverse)`.

**Yang–Gu:** let ``W=G^*G=D+W_{\mathrm{off}}``, where ``D`` is its diagonal.
At each outer step hold ``v=P_b(Gx_k)`` fixed and perform the inner iteration

```math
x^{(m+1)}=P_a\left(G^*v-W_{\mathrm{off}}x^{(m)}\right).
```

At a vanishing argument the previous phase is retained. For ``G^*G=cI`` the
off-diagonal part vanishes and the method reduces to GS. The implementation
computes ``W x`` via forward and adjoint applications, without storing a dense
Gram matrix. The matrix interface uses Julia column-major vectorization.
For a matrix-free operator supply `projector=(forward, adjoint)` and
`normaldiag=diag(G'G)` reshaped like `Ain` (or a scalar for a constant diagonal).
Using the inverse in place of the adjoint is incorrect for a nonunitary system.
This is the unweighted simultaneous inner iteration; poorly conditioned
operators can cause stagnation or oscillation. `inneriter` and `innertol`
control the inner loop separately from the outer measurement error.

**HIO:** with ``p=P_M(x_k)``, accept ``p`` inside support ``S`` and apply
``x_k-\beta p`` outside. With `nonnegative=true`, accepted pixels must also have
``\Re(p)\geq0`` and are replaced by ``\Re(p)``; other pixels receive feedback.
This option models a real nonnegative object; support alone permits complex fields.

**RAAR:** with ``R_S=2P_S-I`` and ``R_M=2P_M-I``, use

```math
x_{k+1}=\frac{\beta}{2}(R_S R_M x_k+x_k)+(1-\beta)P_Mx_k.
```

``P_S`` zeros the exterior of the support and optionally projects its interior
onto the nonnegative real axis. It acts on the reflected field ``R_M x_k``.
At ``\beta=0`` this becomes ``P_M``; at ``\beta=1`` it becomes averaged
alternating reflection. HIO and RAAR retain feedback after the final update.
They do not silently append error-reduction steps.

**SBMIR:** initialize the first plane with ``\sqrt{I_1}`` and the initial phase.
Propagate successively to each plane, replacing amplitude by ``\sqrt{I_j}``.
Then traverse the planes backwards, enforcing each visited amplitude again.
One iteration is the complete forward/backward sweep. For two planes this is
GS with the specified propagation pair. Custom adjacent-plane pairs are passed
as `propagators=[(forward12, backward21), (forward23, backward32), ...]`.

## Fresnel propagation and a runnable example

`fresnelpropagate` implements the paraxial transfer function

```math
H_z(f_1,\ldots,f_d)=\exp\left(-i\pi\lambda z\sum_j f_j^2\right).
```

The constant carrier phase is omitted. Positive distance uses this sign;
negative distance is the inverse. The grid is periodic, fixed between planes,
and may have different spacings along different dimensions. Wavelength,
positions and pixel spacing must share a length unit. Choose adequate sampling
and padding for the optical field. This model assumes paraxial propagation.

```julia
using PhaseRetrieval, Random
Random.seed!(104)
truth = (0.6 .+ rand(8, 9)) .* exp.(0.7im .* randn(8, 9))
z = [0.0, 0.8, 2.1, 3.9, 6.2]
optics = (wavelength=0.08, pixelsize=(0.2, 0.25))
I = [abs2.(fresnelpropagate(truth, d; optics...)) for d in z]
fields = sbmir(I, z; optics..., init=ones(size(truth)), maxiter=250)
phase = angle.(fields[1])
errors = [sum(abs2, abs2.(u) .- measured) / length(u)
          for (u, measured) in zip(fields, I)]
```

These methods solve nonconvex problems: convergence and uniqueness are not
guaranteed. Phase comparisons should account for a global phase offset, and
support-only Fourier retrieval can additionally have translation and conjugate
inversion ambiguities. The tests separately verify update formulas, fixed
points, reduction to simpler algorithms, and recovery from independent initial
estimates; fixed-point tests alone cannot establish recovery performance.

## References

- R. W. Gerchberg and W. O. Saxton, “A practical algorithm for the determination
  of the phase from image and diffraction plane pictures,” Optik **35**, 237–246 (1972).
- G.-Z. Yang et al., [“Gerchberg–Saxton and Yang–Gu algorithms for phase retrieval
  in a nonunitary transform system: a comparison”](https://doi.org/10.1364/AO.33.000209),
  Applied Optics **33**, 209–218 (1994).
- J. R. Fienup, [“Phase retrieval algorithms: a comparison”](https://doi.org/10.1364/AO.21.002758),
  Applied Optics **21**, 2758–2769 (1982).
- D. R. Luke, [“Relaxed averaged alternating reflections for diffraction imaging”](https://doi.org/10.1088/0266-5611/21/1/004),
  Inverse Problems **21**, 37–50 (2005).
- P. Almoro, G. Pedrini and W. Osten, [“Complete wavefront reconstruction using
  sequential intensity measurements of a volume speckle field”](https://doi.org/10.1364/AO.45.008596),
  Applied Optics **45**, 8596–8605 (2006).

## API

```@index
```

```@autodocs
Modules = [PhaseRetrieval]
```
