# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ToZerosBridge{T} <: MOI.Bridges.Constraint.AbstractBridge

`ToZerosBridge` implements the following reformulation:

  * `[x₁, x₂]` in `ComplementsWithSetType{Zeros}` into `[x₁]` in
    [`MOI.Zeros`](@ref)

`ComplementsWithSetType{Zeros}` means that the activity `x₁` is complementary
to a variable `x₂` whose reference set is a single point (it is either free or
fixed). In both cases the complementarity condition forces the activity to be
zero, so the constraint reduces to `x₁ = 0`.

## Source node

`ToZerosBridge` supports:

  * [`MOI.VectorOfVariables`](@ref) in
    [`ComplementsWithSetType{MOI.Zeros}`](@ref)

## Target nodes

`ToZerosBridge` creates:

  * [`MOI.VectorOfVariables`](@ref) in [`MOI.Zeros`](@ref)

"""
struct ToZerosBridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    zeros::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Zeros}
    func::MOI.VectorOfVariables
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToZerosBridge{T}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{MOI.Zeros},
) where {T}
    @assert set.dimension == 2
    x1 = func.variables[1]
    ci = MOI.add_constraint(model, MOI.VectorOfVariables([x1]), MOI.Zeros(1))
    return ToZerosBridge{T}(ci, func)
end

function MOI.supports_constraint(
    ::Type{<:ToZerosBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{MOI.Zeros}},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ToZerosBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{MOI.Zeros}},
) where {T}
    return ToZerosBridge{T}
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:ToZerosBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{ToZerosBridge{T}}) where {T}
    return Tuple{Type,Type}[(MOI.VectorOfVariables, MOI.Zeros)]
end

function MOI.get(
    ::ToZerosBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,MOI.Zeros},
)::Int64
    return 1
end

function MOI.get(
    bridge::ToZerosBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.Zeros},
)
    return [bridge.zeros]
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::ToZerosBridge,
)
    return bridge.func
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, ::ToZerosBridge)
    return ComplementsWithSetType{MOI.Zeros}(2)
end

function MOI.delete(model::MOI.ModelLike, bridge::ToZerosBridge)
    MOI.delete(model, bridge.zeros)
    return
end
