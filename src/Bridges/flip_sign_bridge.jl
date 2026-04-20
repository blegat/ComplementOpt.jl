"""
    FlipSignBridge{T,F,S1,S2} <: MOI.Bridges.Constraint.AbstractBridge

`FlipSignBridge` implements the following reformulations:

  * `[x, y]` in `ComplementsWithSetType{Nonnegatives}` into
    `[-x, y]` in `ComplementsWithSetType{Nonpositives}`
  * `[x, y]` in `ComplementsWithSetType{Nonpositives}` into
    `[-x, y]` in `ComplementsWithSetType{Nonnegatives}`
  * `[x, y]` in `ComplementsWithSetType{GreaterThan{T}}` into
    `[-x, y]` in `ComplementsWithSetType{LessThan{T}}`
  * `[x, y]` in `ComplementsWithSetType{LessThan{T}}` into
    `[-x, y]` in `ComplementsWithSetType{GreaterThan{T}}`

Only the first component (the activity) is negated; the second component
(the slack variable) is unchanged.

## Source node

`FlipSignBridge` supports:

  * `G` in `ComplementsWithSetType{S1}`

## Target nodes

`FlipSignBridge` creates:

  * `F` in `ComplementsWithSetType{S2}`

"""
struct FlipSignBridge{T,F,S1,S2,G} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{F,ComplementsWithSetType{S2}}
    original_func::G
end

function _flip_set_type(::Type{MOI.Nonnegatives})
    return MOI.Nonpositives
end
function _flip_set_type(::Type{MOI.Nonpositives})
    return MOI.Nonnegatives
end
function _flip_set_type(::Type{MOI.GreaterThan{T}}) where {T}
    return MOI.LessThan{T}
end
function _flip_set_type(::Type{MOI.LessThan{T}}) where {T}
    return MOI.GreaterThan{T}
end

const _FlippableSets = Union{MOI.Nonnegatives,MOI.Nonpositives,MOI.GreaterThan,MOI.LessThan}

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{FlipSignBridge{T,F,S1,S2,G}},
    model::MOI.ModelLike,
    func::G,
    set::ComplementsWithSetType{S1},
) where {T,F,S1,S2,G<:MOI.AbstractVectorFunction}
    @assert set.dimension == 2
    scalars = MOIU.scalarize(func)
    # Negate only the first component (the activity)
    neg_x = MOIU.operate(-, T, scalars[1])
    new_func = MOIU.operate(vcat, T, neg_x, scalars[2])
    ci = MOI.add_constraint(model, new_func, ComplementsWithSetType{S2}(2))
    return FlipSignBridge{T,F,S1,S2,G}(ci, func)
end

function MOI.supports_constraint(
    ::Type{<:FlipSignBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {T,S<:_FlippableSets}
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:FlipSignBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S1}},
) where {T,S1<:_FlippableSets}
    H = MOIU.promote_operation(-, T, MOIU.scalar_type(G))
    F = MOIU.promote_operation(vcat, T, H, MOIU.scalar_type(G))
    S2 = _flip_set_type(S1)
    return FlipSignBridge{T,F,S1,S2,G}
end

MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:FlipSignBridge},
) = true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::FlipSignBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.constraint, value)
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:FlipSignBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:FlipSignBridge{T,F,S1,S2}},
) where {T,F,S1,S2}
    return Tuple{Type,Type}[(F, ComplementsWithSetType{S2})]
end

function MOI.get(
    ::FlipSignBridge{T,F,S1,S2},
    ::MOI.NumberOfConstraints{F,ComplementsWithSetType{S2}},
)::Int64 where {T,F,S1,S2}
    return 1
end

function MOI.get(
    bridge::FlipSignBridge{T,F,S1,S2},
    ::MOI.ListOfConstraintIndices{F,ComplementsWithSetType{S2}},
) where {T,F,S1,S2}
    return [bridge.constraint]
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::FlipSignBridge,
)
    return bridge.original_func
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    ::FlipSignBridge{T,F,S1},
) where {T,F,S1}
    return ComplementsWithSetType{S1}(2)
end

function MOI.delete(model::MOI.ModelLike, bridge::FlipSignBridge)
    MOI.delete(model, bridge.constraint)
    return
end
