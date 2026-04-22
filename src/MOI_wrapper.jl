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

"""
    _inner_supports_nlp(model::Optimizer)

Check whether the inner solver supports nonlinear constraints
(`ScalarNonlinearFunction`-in-`LessThan{Float64}`). If true, the nonlinear
relaxation path is used. Otherwise, the SOS1 path is used.
"""
function _inner_supports_nlp(model::Optimizer)
    return MOI.supports_constraint(
        model.model,
        MOI.ScalarNonlinearFunction,
        MOI.LessThan{Float64},
    )
end

# No variable bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractSet}) = false

# Complements and ComplementsWithSetType are bridged
MOI.Bridges.is_bridged(::Optimizer, ::Type{MOI.Complements}) = true
MOI.Bridges.supports_bridging_constrained_variable(::Optimizer, ::Type{MOI.Complements}) =
    true
MOI.Bridges.bridge_type(::Optimizer, ::Type{MOI.Complements}) =
    Bridges.SpecifySetTypeBridge{Float64}

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
) = Bridges.VerticalBridge{Float64,MOI.Complements}

# VectorOfVariables-in-Complements → Bridges.SpecifySetTypeBridge{Float64}
MOI.Bridges.bridge_type(
    ::Optimizer,
    ::Type{<:MOI.VectorOfVariables},
    ::Type{MOI.Complements},
) = Bridges.SpecifySetTypeBridge{Float64}

# ComplementsWithSetType{S} → bridge selection depends on inner solver
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

# --- NLP path: inner solver supports ScalarNonlinearFunction ---
# NonlinearBridge handles all set types directly via relaxation methods.

# --- SOS1 path: inner solver does NOT support ScalarNonlinearFunction ---
# The chain is:
#   ComplementsWithSetType{Interval} → SplitIntervalBridge → {GreaterThan, LessThan}
#   ComplementsWithSetType{LessThan/Nonpositives} → FlipSignBridge → {GreaterThan/Nonnegatives}
#   ComplementsWithSetType{GreaterThan} → ComplementsVectorizeBridge → VAF-in-{Nonnegatives}
#   VAF-in-ComplementsWithSetType{Nonnegatives} → VerticalBridge → VOV-in-{Nonnegatives}
#   VOV-in-ComplementsWithSetType{Nonnegatives} → ToSOS1Bridge → SOS1

function MOI.Bridges.bridge_type(
    b::Optimizer,
    ::Type{<:MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {S}
    if _inner_supports_nlp(b)
        return Bridges.NonlinearBridge{Float64,S}
    end
    return _sos1_bridge_type(MOI.VectorOfVariables, S)
end

function MOI.Bridges.bridge_type(
    b::Optimizer,
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {S}
    if _inner_supports_nlp(b)
        return Bridges.NonlinearBridge{Float64,S}
    end
    return _sos1_bridge_type(F, S)
end

# Dispatch to the appropriate bridge for the SOS1 path.
# The goal is to reach VOV-in-ComplementsWithSetType{Nonnegatives} → ToSOS1Bridge.

# Interval → SplitInterval (split into GreaterThan + LessThan)
_sos1_bridge_type(::Type{MOI.VectorOfVariables}, ::Type{<:MOI.Interval}) =
    Bridges.SplitIntervalBridge{Float64}

# LessThan/Nonpositives → FlipSign (negate activity to get GreaterThan/Nonnegatives)
_sos1_bridge_type(
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.LessThan,MOI.Nonpositives}},
) = Bridges.FlipSignBridge{Float64}

# VOV-in-GreaterThan/EqualTo → Vectorize (shift to Nonneg/Zeros)
_sos1_bridge_type(
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.GreaterThan,MOI.EqualTo}},
) = Bridges.ComplementsVectorizeBridge{Float64}

# VOV-in-Nonnegatives ��� ToSOS1Bridge (final target)
_sos1_bridge_type(::Type{MOI.VectorOfVariables}, ::Type{MOI.Nonnegatives}) =
    Bridges.ToSOS1Bridge{Float64}

# VOV-in-Zeros → ToSOS1Bridge (trivial complementarity)
_sos1_bridge_type(::Type{MOI.VectorOfVariables}, ::Type{MOI.Zeros}) =
    Bridges.ToSOS1Bridge{Float64}

# Any non-VOV function → VerticalBridge (create slacks, then re-enter as VOV)
function _sos1_bridge_type(::Type{F}, ::Type{S}) where {F<:MOI.AbstractVectorFunction,S}
    return Bridges.VerticalBridge{Float64,ComplementsWithSetType{S}}
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
