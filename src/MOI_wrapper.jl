struct Reformulation <: MOI.AbstractOptimizerAttribute end

struct Bridges
end

Base.isempty(::Bridges) = true
MOI.Bridges.Constraint.has_bridges(::Bridges) = false

struct Optimizer{O<:MOI.ModelLike} <: MOI.Bridges.AbstractBridgeOptimizer
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
        )
    end
end

MOI.Bridges.Constraint.bridges(model::Optimizer) = model.constraint_map

# No variable bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractSet}) = false

# No objective bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractFunction}) = false

# We only bridge complements constraints
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractFunction}, ::Type{<:MOI.AbstractSet}) = false
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractVectorFunction}, ::Type{<:MOI.Complements}) = true
MOI.Bridges.supports_bridging_constraint(::Optimizer, ::Type{<:MOI.AbstractVectorFunction}, ::Type{<:MOI.Complements}) = true
MOI.Bridges.bridge_type(::Optimizer, ::Type{<:MOI.AbstractVectorFunction}, ::Type{<:MOI.Complements}) = VerticalBridge
MOI.Bridges.bridge_type(model::Optimizer, ::Type{<:MOI.VectorOfVariables}, ::Type{<:MOI.Complements}) = NonlinearBridge{model.reformulation}

function MOI.Bridges.bridging_cost(b::Optimizer, args...)
    return MOI.Bridges.bridging_cost(MOI.Bridges.bridge_type(b, args...))
end

# We may have a chain of bridges
MOI.Bridges.recursive_model(b::Optimizer) = b
