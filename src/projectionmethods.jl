function _support(support)
    support isa AbstractArray{Bool} || throw(ArgumentError("support must be a Boolean array"))
    isempty(support) && throw(ArgumentError("support must not be empty"))
    Base.require_one_based_indexing(support)
    support
end

function _object_projection(x, support, nonnegative)
    if nonnegative
        return ifelse.(support, complex.(max.(real.(x), 0)), zero(eltype(x)))
    end
    ifelse.(support, x, zero(eltype(x)))
end

function _support_solver(method, Aout, support; beta=0.9, nonnegative=false,
                         projector=(fft, ifft), maxiter=100, errorfunc=msd,
                         tol=eps(), init=nothing)
    _measurement(Aout, "Aout")
    _support(support)
    _iterations(maxiter, tol)
    beta isa Real && isfinite(beta) && 0 <= beta <= 1 ||
        throw(ArgumentError("beta must be in [0, 1]"))
    nonnegative isa Bool || throw(ArgumentError("nonnegative must be Boolean"))
    forward, backward = _pair(projector)
    # Projecting a supplied initial iterate would change the HIO/RAAR recurrence.
    x = _initial_field(zeros(typeof(float(zero(eltype(Aout)))), size(support)), init)
    if init === nothing
        x = _object_projection(x, support, nonnegative)
    end
    y = _propagate(forward, x, Aout)
    for _ in 1:maxiter
        if tol > 0 && errorfunc(abs.(y), Aout) <= tol &&
                      errorfunc(x, _object_projection(x, support, nonnegative)) <= tol
            break
        end
        pm = _propagate(backward, _modulus(Aout, y), support)
        if method === :hio
            if nonnegative
                # Real-object HIO: remove imaginary leakage on accepted pixels.
                accepted = support .& (real.(pm) .>= 0)
                x = ifelse.(accepted, complex.(real.(pm)), x .- beta .* pm)
            else
                x = ifelse.(support, pm, x .- beta .* pm)
            end
        else
            reflected = 2 .* pm .- x
            ps = _object_projection(reflected, support, nonnegative)
            # beta/2 * (R_S R_M + I) + (1-beta) * P_M
            x = beta .* (ps .+ x) .+ (1 - 2beta) .* pm
        end
        y = _propagate(forward, x, Aout)
    end
    x, y
end

"""
    hio(Aout, support; beta=0.9, nonnegative=false, projector=(fft, ifft),
        maxiter=100, errorfunc=msd, tol=eps(), init=nothing)

Fienup hybrid input–output phase retrieval from a measured **amplitude** and
a Boolean object support. With `p = inverse(Aout .* phase(forward(x)))`, accept
`p` inside the support and use `x - beta*p` outside. For `nonnegative=true`,
accept `real(p)` only inside the support where `real(p) >= 0`; use feedback on
all remaining pixels. This option assumes a real, nonnegative object.

Return the raw final iterate and its propagation `(x, forward(x))`. Feedback
outside the support is deliberately retained. To form a support-constrained
estimate use `ifelse.(support, x, zero(eltype(x)))` (and clamp its real part if nonnegative).
`init` is copied without projecting; the default random field is projected onto
the object constraint. The forward/inverse pair must be non-mutating and scaled
unitary. Stop only when both the predicted amplitude MSD and object-projection
MSD are at most `tol`. `errorfunc` replaces MSD in both tests. `tol=0` forces
`maxiter` updates; `maxiter=0` returns the initial iterate. `0 <= beta <= 1`.

Reference: J. R. Fienup, Applied Optics **21**, 2758–2769 (1982),
[doi:10.1364/AO.21.002758](https://doi.org/10.1364/AO.21.002758).
"""
hio(Aout, support; kwargs...) = _support_solver(:hio, Aout, support; kwargs...)

"""
    raar(Aout, support; beta=0.9, nonnegative=false, projector=(fft, ifft),
         maxiter=100, errorfunc=msd, tol=eps(), init=nothing)

Relaxed averaged alternating reflections using
`x_next = beta/2 * (R_S(R_M(x)) + x) + (1-beta)*P_M(x)`,
where `R = 2P-I`, `P_M` imposes the measured amplitude, and `P_S` imposes the
Boolean support and optionally real nonnegativity. In particular, nonnegativity
is applied to the **reflected** field, not to `P_M(x)`.

Return `(x, forward(x))` with the raw iterate, retaining relaxation feedback.
Arguments, stopping, initialization, and support-estimate extraction follow
[`hio`](@ref). `beta=0` is the modulus projection; `beta=1` is averaged
alternating reflection.

Reference: D. R. Luke, Inverse Problems **21**, 37–50 (2005),
[doi:10.1088/0266-5611/21/1/004](https://doi.org/10.1088/0266-5611/21/1/004).
"""
raar(Aout, support; kwargs...) = _support_solver(:raar, Aout, support; kwargs...)
