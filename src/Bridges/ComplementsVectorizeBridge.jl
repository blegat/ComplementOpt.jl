# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ComplementsVectorizeBridge{T,F,S,SV} <: MOI.Bridges.Constraint.AbstractBridge

`ComplementsVectorizeBridge` implements the following reformulations:

  * `[x₁, x₂]` in `ComplementsWithSetType{GreaterThan{T}}` into
    `[x₁, x₂ - lb]` in `ComplementsWithSetType{Nonnegatives}`
  * `[x₁, x₂]` in `ComplementsWithSetType{LessThan{T}}` into
    `[x₁, x₂ - ub]` in `ComplementsWithSetType{Nonpositives}`
  * `[x₁, x₂]` in `ComplementsWithSetType{EqualTo{T}}` into
    `[x₁, x₂ - c]` in `ComplementsWithSetType{Zeros}`

where `T` is the coefficient type.

## Source node

`ComplementsVectorizeBridge` supports:

  * [`MOI.VectorOfVariables`](@ref) in [`ComplementsWithSetType{S}`](@ref)
    where `S` is [`MOI.GreaterThan{T}`](@ref), [`MOI.LessThan{T}`](@ref), or
    [`MOI.EqualTo{T}`](@ref)

## Target nodes

`ComplementsVectorizeBridge` creates:

  * `F` in [`ComplementsWithSetType{SV}`](@ref), where `SV` is
    [`MOI.Nonnegatives`](@ref), [`MOI.Nonpositives`](@ref), or
    [`MOI.Zeros`](@ref) depending on the input set type

"""
mutable struct ComplementsVectorizeBridge{T,F,S,SV} <:
               MOI.Bridges.Constraint.AbstractBridge
    # The shifted constraint is created in `final_touch` (once the bounds of `x2`
    # are set) so it is `nothing` until then.
    constraint::Union{Nothing,MOI.ConstraintIndex{F,ComplementsWithSetType{SV}}}
    func::MOI.VectorOfVariables
    set::ComplementsWithSetType{S}
    reformulation::Union{Nothing,AbstractComplementarityRelaxation}
end

_vector_set_type(::Type{<:MOI.GreaterThan}) = MOI.Nonnegatives

_vector_set_type(::Type{<:MOI.LessThan}) = MOI.Nonpositives

_vector_set_type(::Type{<:MOI.EqualTo}) = MOI.Zeros

function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.GreaterThan},
    x2,
) where {T}
    return MOI.Utilities.get_bounds(model, T, x2)[1]
end

function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.LessThan},
    x2,
) where {T}
    return MOI.Utilities.get_bounds(model, T, x2)[2]
end

function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.EqualTo},
    x2,
) where {T}
    return MOI.Utilities.get_bounds(model, T, x2)[1]
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ComplementsVectorizeBridge{T,F,S,SV}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
) where {T,F,S,SV}
    @assert set.dimension == 2
    # The shift depends on the bounds of `x2`, so it is done in `final_touch`.
    return ComplementsVectorizeBridge{T,F,S,SV}(nothing, func, set, nothing)
end

MOI.Bridges.needs_final_touch(::ComplementsVectorizeBridge) = true

function MOI.Bridges.final_touch(
    bridge::ComplementsVectorizeBridge{T,F,S,SV},
    model::MOI.ModelLike,
) where {T,F,S,SV}
    if bridge.constraint !== nothing
        return
    end
    x1 = bridge.func.variables[1]
    x2 = bridge.func.variables[2]
    c = _set_constant(T, model, bridge.set, x2)
    # Shift only x2: [x1, x2] → [x1, x2 - c]
    shifted = MOI.Utilities.operate(vcat, T, one(T) * x1, one(T) * x2 - c)
    bridge.constraint =
        MOI.add_constraint(model, shifted, ComplementsWithSetType{SV}(2))
    if bridge.reformulation !== nothing
        MOI.set(
            model,
            ComplementarityReformulation(),
            bridge.constraint,
            bridge.reformulation,
        )
    end
    return
end

function MOI.supports_constraint(
    ::Type{<:ComplementsVectorizeBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {S<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ComplementsVectorizeBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}
    G = MOI.Utilities.promote_operation(-, T, MOI.VariableIndex, T)
    F = MOI.Utilities.promote_operation(vcat, T, G, G)
    SV = _vector_set_type(S)
    return ComplementsVectorizeBridge{T,F,S,SV}
end

function MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:ComplementsVectorizeBridge},
)
    return true
end

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::ComplementsVectorizeBridge,
    value::AbstractComplementarityRelaxation,
)
    bridge.reformulation = value
    if bridge.constraint !== nothing
        MOI.set(model, attr, bridge.constraint, value)
    end
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ComplementsVectorizeBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{ComplementsVectorizeBridge{T,F,S,SV}},
) where {T,F,S,SV}
    return Tuple{Type,Type}[(F, ComplementsWithSetType{SV})]
end

function MOI.get(
    bridge::ComplementsVectorizeBridge{T,F,S,SV},
    ::MOI.NumberOfConstraints{F,ComplementsWithSetType{SV}},
)::Int64 where {T,F,S,SV}
    return bridge.constraint === nothing ? 0 : 1
end

function MOI.get(
    bridge::ComplementsVectorizeBridge{T,F,S,SV},
    ::MOI.ListOfConstraintIndices{F,ComplementsWithSetType{SV}},
) where {T,F,S,SV}
    CI = MOI.ConstraintIndex{F,ComplementsWithSetType{SV}}
    return bridge.constraint === nothing ? CI[] : CI[bridge.constraint]
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::ComplementsVectorizeBridge,
)
    return bridge.func
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::ComplementsVectorizeBridge,
)
    return bridge.set
end

function MOI.delete(model::MOI.ModelLike, bridge::ComplementsVectorizeBridge)
    if bridge.constraint !== nothing
        MOI.delete(model, bridge.constraint)
    end
    return
end
