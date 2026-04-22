mutable struct Optimizer{T,O<:MOI.ModelLike} <: MOI.Bridges.AbstractBridgeOptimizer
    model::O # This need to be called `model` by convention of `AbstractBridgeOptimizer`
    reformulation::AbstractComplementarityRelaxation
    constraint_map::MOI.Bridges.Constraint.Map
    con_to_name::Dict{MOI.ConstraintIndex,String}
    name_to_con::Union{Dict{String,MOI.ConstraintIndex},Nothing}
    function Optimizer{T}(model::MOI.ModelLike) where {T}
        return new{T,typeof(model)}(
            model,
            ScholtesRelaxation(zero(T)),
            MOI.Bridges.Constraint.Map(),
            Dict{MOI.ConstraintIndex,String}(),
            nothing,
        )
    end
end

Optimizer(model::MOI.ModelLike) = Optimizer{Float64}(model)

MOI.Bridges.Constraint.bridges(model::Optimizer) = model.constraint_map

"""
    _inner_supports_nlp(model::Optimizer{T}) where {T}

Check whether the inner solver supports nonlinear constraints
(`ScalarNonlinearFunction`-in-`LessThan{T}`). If true, the nonlinear
relaxation path is used. Otherwise, the SOS1 path is used.
"""
function _inner_supports_nlp(model::Optimizer{T}) where {T}
    return MOI.supports_constraint(
        model.model,
        MOI.ScalarNonlinearFunction,
        MOI.LessThan{T},
    )
end

# No variable bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractSet}) = false

# Complements and ComplementsWithSetType are bridged
MOI.Bridges.is_bridged(::Optimizer, ::Type{MOI.Complements}) = true
MOI.Bridges.supports_bridging_constrained_variable(::Optimizer, ::Type{MOI.Complements}) =
    true
MOI.Bridges.bridge_type(::Optimizer{T}, ::Type{MOI.Complements}) where {T} =
    Bridges.SpecifySetTypeBridge{T}

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
    ::Optimizer{T},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) where {T} = Bridges.VerticalBridge{T,MOI.Complements}

# VectorOfVariables-in-Complements → Bridges.SpecifySetTypeBridge{T}
MOI.Bridges.bridge_type(
    ::Optimizer{T},
    ::Type{<:MOI.VectorOfVariables},
    ::Type{MOI.Complements},
) where {T} = Bridges.SpecifySetTypeBridge{T}

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
    b::Optimizer{T},
    ::Type{<:MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    if _inner_supports_nlp(b)
        return Bridges.NonlinearBridge{T,S}
    end
    return _sos1_bridge_type(T, MOI.VectorOfVariables, S)
end

function MOI.Bridges.bridge_type(
    b::Optimizer{T},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    if _inner_supports_nlp(b)
        return Bridges.NonlinearBridge{T,S}
    end
    return _sos1_bridge_type(T, F, S)
end

# Dispatch to the appropriate bridge for the SOS1 path.
# The goal is to reach VOV-in-ComplementsWithSetType{Nonnegatives} → ToSOS1Bridge.

# Interval → SplitInterval (split into GreaterThan + LessThan)
_sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:MOI.Interval},
) where {T} = Bridges.SplitIntervalBridge{T}

# LessThan/Nonpositives → FlipSign (negate activity to get GreaterThan/Nonnegatives)
_sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.LessThan,MOI.Nonpositives}},
) where {T} = Bridges.FlipSignBridge{T}

# VOV-in-GreaterThan/EqualTo → Vectorize (shift to Nonneg/Zeros)
_sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.GreaterThan,MOI.EqualTo}},
) where {T} = Bridges.ComplementsVectorizeBridge{T}

# VOV-in-Nonnegatives → ToSOS1Bridge (final target)
_sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Nonnegatives},
) where {T} = Bridges.ToSOS1Bridge{T}

# VOV-in-Zeros → ToSOS1Bridge (trivial complementarity)
_sos1_bridge_type(::Type{T}, ::Type{MOI.VectorOfVariables}, ::Type{MOI.Zeros}) where {T} =
    Bridges.ToSOS1Bridge{T}

# Any non-VOV function → VerticalBridge (create slacks, then re-enter as VOV)
function _sos1_bridge_type(
    ::Type{T},
    ::Type{F},
    ::Type{S},
) where {T,F<:MOI.AbstractVectorFunction,S}
    return Bridges.VerticalBridge{T,ComplementsWithSetType{S}}
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
