mutable struct Optimizer{O<:MOI.ModelLike} <: MOI.Bridges.AbstractBridgeOptimizer
    model::O # This need to be called `model` by convention of `AbstractBridgeOptimizer`
    reformulation::AbstractComplementarityRelaxation
    constraint_map::MOI.Bridges.Constraint.Map
    con_to_name::Dict{MOI.ConstraintIndex,String}
    name_to_con::Union{Dict{String,MOI.ConstraintIndex},Nothing}
    function Optimizer(model::MOI.ModelLike)
        return new{typeof(model)}(
            model,
            ScholtesRelaxation(0.0),
            MOI.Bridges.Constraint.Map(),
            Dict{MOI.ConstraintIndex,String}(),
            nothing,
        )
    end
end

MOI.Bridges.Constraint.bridges(model::Optimizer) = model.constraint_map

# No variable bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractSet}) = false

# Complements and ComplementsWithSetType are bridged
MOI.Bridges.is_bridged(::Optimizer, ::Type{MOI.Complements}) = true
MOI.Bridges.supports_bridging_constrained_variable(::Optimizer, ::Type{MOI.Complements}) =
    true
MOI.Bridges.bridge_type(::Optimizer, ::Type{MOI.Complements}) = SpecifySetTypeBridge{Float64}

MOI.Bridges.is_bridged(::Optimizer, ::Type{<:ComplementsWithSetType}) = true
MOI.Bridges.supports_bridging_constrained_variable(
    ::Optimizer,
    ::Type{<:ComplementsWithSetType},
) = true

# No objective bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractFunction}) = false

# We only bridge Complements and ComplementsWithSetType constraints
MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractFunction},
    ::Type{<:MOI.AbstractSet},
) = false

# Expression-based Complements → Bridges.VerticalBridge
MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) = true
MOI.Bridges.supports_bridging_constraint(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) = true
MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) = Bridges.VerticalBridge{MOI.Complements}

# VectorOfVariables-in-Complements → Bridges.SpecifySetTypeBridge{Float64}
MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.VectorOfVariables},
    ::Type{MOI.Complements},
) = Bridges.SpecifySetTypeBridge{Float64}

# ComplementsWithSetType{S} → Bridges.NonlinearBridge{S} for all S
MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
) = true
MOI.Bridges.supports_bridging_constraint(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
) = true

function MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {S}
    return Bridges.NonlinearBridge{S}
end

function MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {S}
    return Bridges.NonlinearBridge{S}
end

function MOI.Bridges.bridging_cost(b::Optimizer, args...)
    return MOI.Bridges.bridging_cost(MOI.Bridges.bridge_type(b, args...))
end

# We may have a chain of bridges
MOI.Bridges.recursive_model(b::Optimizer) = b

MOI.supports(::Optimizer, ::DefaultComplementarityReformulation) = true

function MOI.set(
    model::Optimizer,
    ::DefaultComplementarityReformulation,
    reformulation::AbstractComplementarityRelaxation,
)
    model.reformulation = reformulation
    return
end

MOI.Utilities.map_indices(::Function, relax::AbstractComplementarityRelaxation) = relax

_additional_arguments(::Optimizer, ::Type) = tuple()

function _additional_arguments(model::Optimizer, ::Type{<:Bridges.NonlinearBridge})
    # Create a 1-element tuple since it is splatted in `add_bridged_constraint`
    return (model.reformulation,)
end

# TODO it would be nice if MOI was defining this `MOI.Bridges.additional_arguments` function and
#      already had this implementation of `add_bridged_constraint` so that I don't have to reimplement it
function MOI.Bridges.add_bridged_constraint(b::Optimizer, BridgeType, f, s)
    bridge = MOI.Bridges.Constraint.Constraint.bridge_constraint(
        BridgeType,
        MOI.Bridges.recursive_model(b),
        f,
        s,
        _additional_arguments(b, BridgeType)...,
    )
    # The rest is copy-pasted from the default implementation of `add_bridged_constraint` in MOI
    ci = MOI.Bridges.Constraint.add_key_for_bridge(
        MOI.Bridges.Constraint.bridges(b)::MOI.Bridges.Constraint.Map,
        bridge,
        f,
        s,
        !Base.Fix1(MOI.is_valid, MOI.Bridges.Variable.bridges(b)),
    )
    MOI.Bridges.Variable.register_context(MOI.Bridges.Variable.bridges(b), ci)
    return ci
end
