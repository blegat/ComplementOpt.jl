# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

mutable struct Optimizer{T,O<:MOI.ModelLike} <:
               MOI.Bridges.AbstractBridgeOptimizer
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

MOI.Bridges.supports_constraint_bridges(::Optimizer) = true
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

function MOI.Bridges.supports_bridging_constrained_variable(
    ::Optimizer,
    ::Type{MOI.Complements},
)
    return true
end

function MOI.Bridges.bridge_type(
    ::Optimizer{T},
    ::Type{MOI.Complements},
) where {T}
    return Bridges.SpecifySetTypeBridge{T}
end

MOI.Bridges.is_bridged(::Optimizer, ::Type{<:ComplementsWithSetType}) = true

function MOI.Bridges.supports_bridging_constrained_variable(
    ::Optimizer,
    ::Type{<:ComplementsWithSetType},
)
    return true
end

# No objective bridge
MOI.Bridges.is_bridged(::Optimizer, ::Type{<:MOI.AbstractFunction}) = false

# We only bridge Complements and ComplementsWithSetType constraints
function MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractFunction},
    ::Type{<:MOI.AbstractSet},
)
    return false
end

# Expression-based Complements → Bridges.VerticalBridge
function MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
)
    return true
end

function MOI.Bridges.supports_bridging_constraint(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
)
    return true
end

function MOI.Bridges.bridge_type(
    ::Optimizer{T},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) where {T}
    return Bridges.VerticalBridge{T,MOI.Complements}
end

# VectorOfVariables-in-Complements → Bridges.SpecifySetTypeBridge{T}
function MOI.Bridges.bridge_type(
    ::Optimizer{T},
    ::Type{<:MOI.VectorOfVariables},
    ::Type{MOI.Complements},
) where {T}
    return Bridges.SpecifySetTypeBridge{T}
end

# ComplementsWithSetType{S} → bridge selection depends on inner solver
function MOI.Bridges.is_bridged(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

function MOI.Bridges.supports_bridging_constraint(
    ::Optimizer,
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

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
function _sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:MOI.Interval},
) where {T}
    return Bridges.SplitIntervalBridge{T}
end

# LessThan/Nonpositives → FlipSign (negate activity to get GreaterThan/Nonnegatives)
function _sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.LessThan,MOI.Nonpositives}},
) where {T}
    return Bridges.FlipSignBridge{T}
end

# VOV-in-GreaterThan/EqualTo → Vectorize (shift to Nonneg/Zeros)
function _sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.GreaterThan,MOI.EqualTo}},
) where {T}
    return Bridges.ComplementsVectorizeBridge{T}
end

# VOV-in-Nonnegatives → ToSOS1Bridge (final target)
function _sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Nonnegatives},
) where {T}
    return Bridges.ToSOS1Bridge{T}
end

# VOV-in-Zeros → ToZerosBridge (trivial complementarity: activity is zero)
function _sos1_bridge_type(
    ::Type{T},
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Zeros},
) where {T}
    return Bridges.ToZerosBridge{T}
end

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

function MOI.Utilities.map_indices(
    ::Function,
    relax::AbstractComplementarityRelaxation,
)
    return relax
end

_additional_arguments(::Optimizer, ::Type) = tuple()

function _additional_arguments(
    model::Optimizer,
    ::Type{<:Bridges.NonlinearBridge},
)
    # Create a 1-element tuple since it is splatted in `add_bridged_constraint`
    return (model.reformulation,)
end

# TODO(blegat): it would be nice if MOI was defining this
# `MOI.Bridges.additional_arguments` function and already had this
# implementation of `add_bridged_constraint` so that I don't have to reimplement
# it
function MOI.Bridges.add_bridged_constraint(b::Optimizer, BridgeType, f, s)
    bridge = MOI.Bridges.Constraint.Constraint.bridge_constraint(
        BridgeType,
        MOI.Bridges.recursive_model(b),
        f,
        s,
        _additional_arguments(b, BridgeType)...,
    )
    # The rest is copy-pasted from the default implementation of
    # `add_bridged_constraint` in MOI
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
