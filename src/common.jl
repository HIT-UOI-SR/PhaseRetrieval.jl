# Public solvers and user-supplied propagation functions are non-mutating.
function _measurement(a, name)
    a isa AbstractArray{<:Real} || throw(ArgumentError("$name must be a real array"))
    isempty(a) && throw(ArgumentError("$name must not be empty"))
    Base.require_one_based_indexing(a)
    all(x -> isfinite(x) && x >= 0, a) ||
        throw(ArgumentError("$name must contain finite, nonnegative values"))
    a
end

function _iterations(maxiter, tol)
    maxiter isa Integer && !(maxiter isa Bool) && maxiter >= 0 ||
        throw(ArgumentError("maxiter must be a nonnegative integer"))
    tol isa Real && isfinite(tol) && tol >= 0 ||
        throw(ArgumentError("tol must be finite and nonnegative"))
end

_same_size(a, b, name) =
    size(a) == size(b) || throw(DimensionMismatch("$name has size $(size(a)); expected $(size(b))"))

function _initial_field(a, init)
    T = typeof(real(float(zero(eltype(a)))))
    if init === nothing
        return randn(Complex{T}, size(a))
    end
    init isa AbstractArray{<:Number} || throw(ArgumentError("init must be a numeric array"))
    Base.require_one_based_indexing(init)
    _same_size(init, a, "init")
    all(isfinite, init) || throw(ArgumentError("init must contain finite values"))
    R = promote_type(T, typeof(real(float(zero(eltype(init))))))
    Complex{R}.(init)
end

# The support fixes the object shape. Light fields also need physical metadata.
_support_initial_field(a, support, init) =
    _initial_field(zeros(typeof(real(float(zero(eltype(a))))), size(support)), init)

# At zero choose phase zero, so a positive target amplitude is not lost.
_unitphase(z) = iszero(z) ? one(complex(z)) : z / abs(z)
_modulus(a, z) = a .* _unitphase.(z)

function _propagate(f, x, target)
    y = f(x)
    y isa AbstractArray{<:Number} || throw(ArgumentError("propagator must return a numeric array"))
    Base.require_one_based_indexing(y)
    _same_size(y, target, "propagator output")
    all(isfinite, y) || throw(ArgumentError("propagator returned nonfinite values"))
    y
end

function _pair(projector)
    length(projector) == 2 || throw(ArgumentError("projector must contain (forward, backward)"))
    projector[1], projector[2]
end

# OpticalPropagation broadcasts retain a ComplexF64 light-field container even
# for abs/abs2. Validate the values without discarding the physical grid.
function _measurement(a::MonoLightField2D, name)
    all(isreal, a) || throw(ArgumentError("$name must contain real measurements"))
    _measurement(real.(parent(a)), name)
    a
end

_initial_field(a::MonoLightField2D, ::Nothing) =
    MonoLightField2D(a; data=invoke(_initial_field, Tuple{Any,Any}, a, nothing))

function _support_initial_field(a::MonoLightField2D, support, init)
    if init === nothing
        # A different object grid cannot be inferred from a detector grid.
        _same_size(support, a, "support for default light-field initialization")
        return _initial_field(a, nothing)
    end
    invoke(_support_initial_field, Tuple{Any,Any,Any}, a, support, init)
end
