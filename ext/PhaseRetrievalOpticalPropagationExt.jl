module PhaseRetrievalOpticalPropagationExt

using PhaseRetrieval
using OpticalPropagation: MonoLightField2D

# OpticalPropagation broadcasts retain a ComplexF64 light-field container even
# for abs/abs2. Validate the values without discarding the physical grid.
function PhaseRetrieval._measurement(a::MonoLightField2D, name)
    all(isreal, a) || throw(ArgumentError("$name must contain real measurements"))
    PhaseRetrieval._measurement(real.(parent(a)), name)
    a
end

function PhaseRetrieval._initial_field(a::MonoLightField2D, init)
    field = invoke(PhaseRetrieval._initial_field, Tuple{Any,Any}, a, init)
    init === nothing ? MonoLightField2D(a; data=field) : field
end

function PhaseRetrieval._support_initial_field(a::MonoLightField2D, support, init)
    if init === nothing
        # A different object grid cannot be inferred from a detector grid.
        PhaseRetrieval._same_size(support, a, "support for default light-field initialization")
        return PhaseRetrieval._initial_field(a, nothing)
    end
    invoke(PhaseRetrieval._support_initial_field, Tuple{Any,Any,Any}, a, support, init)
end

end
