struct Reformulation <: MOI.AbstractConstraintAttribute end

mutable struct Optimizer{O<:MOI.ModelLike} <: MOI.Bridges.AbstractBridgeOptimizer
    model::O # This need to be called `model` by convention of `AbstractBridgeOptimizer`
    reformulation::AbstractComplementarityRelaxation
    constraint_reformulations::Dict{MOI.ConstraintIndex,AbstractComplementarityRelaxation}
    constraint_map::MOI.Bridges.Constraint.Map
    con_to_name::Dict{MOI.ConstraintIndex,String}
    name_to_con::Union{Dict{String,MOI.ConstraintIndex},Nothing}
    function Optimizer(model::MOI.ModelLike)
        return new{typeof(model)}(
            model,
            ScholtesRelaxation(0.0),
            Dict{MOI.ConstraintIndex,AbstractComplementarityRelaxation}(),
            MOI.Bridges.Constraint.Map(),
            Dict{MOI.ConstraintIndex,String}(),
            nothing,
        )
    end
end

MOI.Bridges.Constraint.bridges(model::Optimizer) = model.constraint_map

# No variable bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractSet}) = false

# No objective bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractFunction}) = false

# We only bridge complements constraints
MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractFunction},
    ::Type{<:MOI.AbstractSet},
) = false
MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:MOI.Complements},
) = true
MOI.Bridges.supports_bridging_constraint(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:MOI.Complements},
) = true
MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:MOI.Complements},
) = VerticalBridge

# It's a bit unfortunate that we're passing the reformulation as type parameter.
# This means that if the user change the **value** of the tolerance in the reformulation
# then it will trigger a recompilation of the bridge.
MOI.Bridges.bridge_type(
    model::Optimizer,
    ::Type{<:MOI.VectorOfVariables},
    ::Type{<:MOI.Complements},
) = NonlinearBridge

function MOI.Bridges.bridging_cost(b::Optimizer, args...)
    return MOI.Bridges.bridging_cost(MOI.Bridges.bridge_type(b, args...))
end

# We may have a chain of bridges
MOI.Bridges.recursive_model(b::Optimizer) = b

# Relaxation attribute
struct RelaxationMethod <: MOI.AbstractOptimizerAttribute end

MOI.supports(::Optimizer, ::RelaxationMethod) = true

function MOI.set(
    model::Optimizer,
    ::RelaxationMethod,
    reformulation::AbstractComplementarityRelaxation,
)
    model.reformulation = reformulation
    return
end

MOI.Utilities.map_indices(::Function, relax::AbstractComplementarityRelaxation) = relax

function MOI.supports(
    ::Optimizer,
    ::Reformulation,
    ::Type{<:MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:MOI.Complements}},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::Reformulation,
    ci::MOI.ConstraintIndex,
    value::AbstractComplementarityRelaxation,
)
    model.constraint_reformulations[ci] = value
    return
end

function MOI.get(
    model::Optimizer,
    ::Reformulation,
    ci::MOI.ConstraintIndex,
)
    return get(model.constraint_reformulations, ci, model.reformulation)
end

# Forward the Reformulation attribute through JuMP's LazyBridgeOptimizer
function MOI.supports(
    ::MOI.Bridges.LazyBridgeOptimizer{<:Optimizer},
    ::Reformulation,
    ::Type{<:MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,<:MOI.Complements}},
)
    return true
end

function MOI.set(
    b::MOI.Bridges.LazyBridgeOptimizer{<:Optimizer},
    attr::Reformulation,
    ci::MOI.ConstraintIndex,
    value::AbstractComplementarityRelaxation,
)
    return MOI.set(b.model, attr, ci, value)
end

function MOI.get(
    b::MOI.Bridges.LazyBridgeOptimizer{<:Optimizer},
    attr::Reformulation,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(b.model, attr, ci)
end

_additional_arguments(::Optimizer, ::Type) = tuple()

function _additional_arguments(model::Optimizer, ::Type{NonlinearBridge})
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

function _rebridge!(model::Optimizer, bridge::NonlinearBridge, reformulation::AbstractComplementarityRelaxation)
    # Delete old constraints from the inner model
    for ci in bridge.constraints
        MOI.delete(model.model, ci)
    end
    # Re-create with the new reformulation
    new_constraints = reformulate_as_nonlinear_program!(
        model.model,
        reformulation,
        bridge.func,
        bridge.set,
    )
    empty!(bridge.constraints)
    append!(bridge.constraints, new_constraints)
    return
end

function _get_nonlinear_bridge(model::Optimizer, ci::MOI.ConstraintIndex)
    bridge = MOI.Bridges.bridge(model, ci)
    if bridge isa VerticalBridge
        return MOI.Bridges.bridge(model, bridge.constraint)
    end
    return bridge
end

function MOI.optimize!(model::Optimizer)
    for (ci, relax) in model.constraint_reformulations
        bridge = _get_nonlinear_bridge(model, ci)
        _rebridge!(model, bridge, relax)
    end
    return MOI.optimize!(model.model)
end
