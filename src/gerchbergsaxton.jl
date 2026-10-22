"""
    gerchbergsaxton(Ain, Aout; maxiter=100, projector=(fft, ifft),
                   errorfunc=msd, tol=eps(), init=nothing)

Recover a complex field from two **amplitude** measurements using Gerchberg–Saxton.
For intensity measurements, pass their square roots. `projector` contains a
non-mutating forward transform and its inverse (normally a scaled unitary pair).
The default FFT uses unshifted frequency ordering and Julia's FFT normalization.

Return `(field, propagated)`, where `abs.(field) ≈ Ain` and
`propagated = projector[1](field)`. The latter is a prediction, not an amplitude
clamped array: its residual against `Aout` measures consistency. Stop when
`errorfunc(abs.(propagated), Aout) <= tol`, or after `maxiter` updates.
`tol=0` disables early stopping. `maxiter=0` returns the amplitude-constrained
initial estimate and its propagation. `init` is a field estimate (only
its phase is used), defaults to random phase, and is never modified.

## References

1. R. W. Gerchberg and W. O. Saxton, "A practical algorithm for the determination of the phase from image and diffraction plane pictures," Optik 35, 237 (1972)
"""
function gerchbergsaxton(Ain, Aout; maxiter=100, projector=(fft, ifft),
                        errorfunc=msd, tol=eps(), init=nothing)
    _measurement(Ain, "Ain")
    _measurement(Aout, "Aout")
    _iterations(maxiter, tol)
    forward, backward = _pair(projector)
    x = _modulus(Ain, _initial_field(Ain, init))
    y = _propagate(forward, x, Aout)
    for _ in 1:maxiter
        tol > 0 && errorfunc(abs.(y), Aout) <= tol && break
        x = _modulus(Ain, _propagate(backward, _modulus(Aout, y), Ain))
        y = _propagate(forward, x, Aout)
    end
    x, y
end
