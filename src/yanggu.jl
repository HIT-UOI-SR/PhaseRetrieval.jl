"""
    yanggu(Ain, Aout, G; kwargs...)
    yanggu(Ain, Aout; projector=nothing, normaldiag=nothing,
           maxiter=100, inneriter=10, innertol=eps(), tol=eps(),
           errorfunc=msd, init=nothing)

Yang–Gu phase retrieval from two **amplitude** measurements for a general linear
transform. A matrix `G` acts on `vec(field)` in Julia column-major order and may
be rectangular or singular. Alternatively, supply `projector=(forward, adjoint)`
and `normaldiag = diag(G'G)` as a scalar or an array shaped like `Ain`.
The second operator is the **Hermitian adjoint**, not the inverse. With neither
keyword supplied, use `(fft, bfft)` and `normaldiag=length(Ain)`.

For each outer step set `v = Aout .* phase(G*u)`. Holding `v` fixed, perform up
to `inneriter` simultaneous inner updates
`u = Ain .* phase(G'*v - (G'*G - Diagonal(normaldiag))*u)`.
At a zero update argument retain the previous phase. This is the unweighted
Yang–Gu inner iteration; convergence is not guaranteed for every nonunitary
operator. For a scaled unitary transform it reduces to Gerchberg–Saxton.
`innertol` bounds the mean squared complex change within the inner loop;
zero disables that stopping test. `inneriter` must be positive.

Return `(field, forward(field))`. Initialization, outer stopping, `tol=0`, and
`maxiter=0` follow [`gerchbergsaxton`](@ref). Inputs are never modified.

Reference: G.-Z. Yang et al., Applied Optics **33**, 209–218 (1994),
[doi:10.1364/AO.33.000209](https://doi.org/10.1364/AO.33.000209).
"""
function yanggu(Ain, Aout; projector=nothing, normaldiag=nothing,
                maxiter=100, inneriter=10, innertol=eps(), tol=eps(),
                errorfunc=msd, init=nothing)
    _measurement(Ain, "Ain")
    _measurement(Aout, "Aout")
    _iterations(maxiter, tol)
    _iterations(inneriter, innertol)
    inneriter > 0 || throw(ArgumentError("inneriter must be positive"))
    if projector === nothing
        projector = (fft, bfft)
        normaldiag === nothing && (normaldiag = length(Ain))
    elseif normaldiag === nothing
        throw(ArgumentError("normaldiag is required with a custom forward/adjoint pair"))
    end
    if normaldiag isa Real
        isfinite(normaldiag) && normaldiag >= 0 ||
            throw(ArgumentError("normaldiag must be finite and nonnegative"))
    else
        _measurement(normaldiag, "normaldiag")
        _same_size(normaldiag, Ain, "normaldiag")
    end
    forward, backward = _pair(projector)
    x = _modulus(Ain, _initial_field(Ain, init))
    y = _propagate(forward, x, Aout)
    for _ in 1:maxiter
        tol > 0 && errorfunc(abs.(y), Aout) <= tol && break
        b = _propagate(backward, _modulus(Aout, y), Ain)
        for _ in 1:inneriter
            normal = _propagate(backward, _propagate(forward, x, Aout), Ain)
            direction = b .- normal .+ normaldiag .* x
            # A vanishing argument leaves the phase undetermined (e.g. a zero column).
            direction = ifelse.(iszero.(direction), x, direction)
            next = _modulus(Ain, direction)
            change = msd(next, x)
            x = next
            innertol > 0 && change <= innertol && break
        end
        y = _propagate(forward, x, Aout)
    end
    x, y
end

function yanggu(Ain, Aout, G::AbstractMatrix; kwargs...)
    _measurement(Ain, "Ain")
    _measurement(Aout, "Aout")
    size(G) == (length(Aout), length(Ain)) ||
        throw(DimensionMismatch("G must have size (length(Aout), length(Ain))"))
    all(isfinite, G) || throw(ArgumentError("G must contain finite values"))
    forward = x -> reshape(G * vec(x), size(Aout))
    backward = y -> reshape(G' * vec(y), size(Ain))
    diagonal = reshape(vec(sum(abs2, G; dims=1)), size(Ain))
    yanggu(Ain, Aout; projector=(forward, backward), normaldiag=diagonal, kwargs...)
end
