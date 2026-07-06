# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    SplitIntervalBridge{T,G} <: MOI.Bridges.Constraint.AbstractBridge

`SplitIntervalBridge` implements the following reformulation:

  * `[x, y]` in `ComplementsWithSetType{Interval{T}}` into
    `[xp, y]` in `ComplementsWithSetType{GreaterThan{T}}` and
    `[xn, y]` in `ComplementsWithSetType{LessThan{T}}`

with the equality constraint `x == xp + xn`, where `xp` and `xn` are new
variables representing the positive and negative parts of the activity.

The input function can be any `MOI.AbstractVectorFunction` (the first
component `x` may be affine or quadratic); only the second component `y`
must be a variable.

Since the lower bound of `y` is strictly smaller than its upper bound,
`xp` and `xn` cannot both be nonzero at a feasible point, so
`xp = max(x, 0)` and `xn = min(x, 0)`. In
[`MOI.Bridges.final_touch`](@ref), if the activity `x` has a finite upper
(resp. lower) bound `ub` (resp. `lb`), the bound `xp <= max(ub, 0)`
(resp. `xn >= min(lb, 0)`) is added. This is done in
[`MOI.Bridges.final_touch`](@ref) and not in
`MOI.Bridges.Constraint.bridge_constraint` because the bounds of `x` may
be set after the constraint is bridged. These bounds are needed by
[`MOI.Bridges.Constraint.SOS1ToMILPBridge`](@ref) in case the inner solver
does not support [`MOI.SOS1`](@ref) constraints.

## Source node

`SplitIntervalBridge` supports:

  * [`MOI.AbstractVectorFunction`](@ref) in
    [`ComplementsWithSetType{MOI.Interval{T}}`](@ref)

## Target nodes

`SplitIntervalBridge` creates:

  * [`MOI.VectorOfVariables`](@ref) in
    [`ComplementsWithSetType{MOI.GreaterThan{T}}`](@ref)
  * [`MOI.VectorOfVariables`](@ref) in
    [`ComplementsWithSetType{MOI.LessThan{T}}`](@ref)
  * `G` in [`MOI.EqualTo{T}`](@ref) (the splitting equality)
  * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{T}`](@ref) and
    [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{T}`](@ref)
    (the bounds on `xp` and `xn`, added in
    [`MOI.Bridges.final_touch`](@ref) if the activity is bounded)

where `G` is the scalar function type of the first component.

"""
mutable struct SplitIntervalBridge{
    T,
    G<:MOI.AbstractScalarFunction,
    F<:MOI.AbstractVectorFunction,
} <: MOI.Bridges.Constraint.AbstractBridge
    lower::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    }
    upper::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    }
    equality::MOI.ConstraintIndex{G,MOI.EqualTo{T}}
    xp::MOI.VariableIndex
    xn::MOI.VariableIndex
    func::F
    set::ComplementsWithSetType{MOI.Interval{T}}
    xp_upper::Union{
        Nothing,
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{T}},
    }
    xn_lower::Union{
        Nothing,
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{T}},
    }
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SplitIntervalBridge{T,G,F}},
    model::MOI.ModelLike,
    func::F,
    set::ComplementsWithSetType{MOI.Interval{T}},
) where {T,G,F<:MOI.AbstractVectorFunction}
    @assert set.dimension == 2
    scalars = MOI.Utilities.scalarize(func)
    x_func = scalars[1]  # activity (may be an expression)
    y = scalars[2]       # slack (must be a variable)
    # y must be a single variable
    y_var = if y isa MOI.VariableIndex
        y
    else
        # Extract the variable from a ScalarAffineFunction wrapping a single variable
        @assert length(y.terms) == 1 &&
                isone(y.terms[1].coefficient) &&
                iszero(y.constant)
        y.terms[1].variable
    end
    # Create xp >= 0 and xn <= 0
    xp, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(zero(T)))
    xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(zero(T)))
    # x == xp + xn
    eq_func = MOI.Utilities.operate(-, T, x_func, xp)
    eq_func = MOI.Utilities.operate!(-, T, eq_func, xn)
    equality = MOI.Utilities.normalize_and_add_constraint(
        model,
        eq_func,
        MOI.EqualTo(zero(T)),
    )
    # [xp, y] in ComplementsWithSetType{GreaterThan{T}}
    lower = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([xp, y_var]),
        ComplementsWithSetType{MOI.GreaterThan{T}}(2),
    )
    # [xn, y] in ComplementsWithSetType{LessThan{T}}
    upper = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([xn, y_var]),
        ComplementsWithSetType{MOI.LessThan{T}}(2),
    )
    return SplitIntervalBridge{T,G,typeof(func)}(
        lower,
        upper,
        equality,
        xp,
        xn,
        func,
        set,
        nothing,
        nothing,
    )
end

function MOI.supports_constraint(
    ::Type{<:SplitIntervalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{MOI.Interval{T}}},
) where {T}
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SplitIntervalBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{MOI.Interval{T}}},
) where {T}
    G = MOI.Utilities.scalar_type(F)
    # After `operate(-, T, G, VariableIndex)`, the type may promote
    H = MOI.Utilities.promote_operation(-, T, G, MOI.VariableIndex)
    return SplitIntervalBridge{T,H,F}
end

function MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:SplitIntervalBridge},
)
    return true
end

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::SplitIntervalBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.lower, value)
    MOI.set(model, attr, bridge.upper, value)
    return
end

function _activity_bounds(model, ::Type{T}, x::MOI.VariableIndex) where {T}
    return MOI.Utilities.get_bounds(model, T, x)
end

function _activity_bounds(
    model,
    ::Type{T},
    func::MOI.ScalarAffineFunction{T},
) where {T}
    f = MOI.Utilities.canonical(func)
    lb = ub = f.constant
    for term in f.terms
        l, u = MOI.Utilities.get_bounds(model, T, term.variable)
        if term.coefficient >= 0
            lb += term.coefficient * l
            ub += term.coefficient * u
        else
            lb += term.coefficient * u
            ub += term.coefficient * l
        end
    end
    return lb, ub
end

function _activity_bounds(
    model,
    ::Type{T},
    ::MOI.AbstractScalarFunction,
) where {T}
    return typemin(T), typemax(T)
end

MOI.Bridges.needs_final_touch(::SplitIntervalBridge) = true

function MOI.Bridges.final_touch(
    bridge::SplitIntervalBridge{T},
    model::MOI.ModelLike,
) where {T}
    x = first(MOI.Utilities.eachscalar(bridge.func))
    lb, ub = _activity_bounds(model, T, x)
    if isfinite(ub)
        set = MOI.LessThan(max(ub, zero(T)))
        if bridge.xp_upper === nothing
            bridge.xp_upper = MOI.add_constraint(model, bridge.xp, set)
        elseif MOI.get(model, MOI.ConstraintSet(), bridge.xp_upper) != set
            MOI.set(model, MOI.ConstraintSet(), bridge.xp_upper, set)
        end
    elseif bridge.xp_upper !== nothing
        MOI.delete(model, bridge.xp_upper)
        bridge.xp_upper = nothing
    end
    if isfinite(lb)
        set = MOI.GreaterThan(min(lb, zero(T)))
        if bridge.xn_lower === nothing
            bridge.xn_lower = MOI.add_constraint(model, bridge.xn, set)
        elseif MOI.get(model, MOI.ConstraintSet(), bridge.xn_lower) != set
            MOI.set(model, MOI.ConstraintSet(), bridge.xn_lower, set)
        end
    elseif bridge.xn_lower !== nothing
        MOI.delete(model, bridge.xn_lower)
        bridge.xn_lower = nothing
    end
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:SplitIntervalBridge{T}},
) where {T}
    return Tuple{Type}[(MOI.GreaterThan{T},), (MOI.LessThan{T},)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:SplitIntervalBridge{T,G}},
) where {T,G}
    return Tuple{Type,Type}[
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.GreaterThan{T}}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.LessThan{T}}),
        (G, MOI.EqualTo{T}),
        (MOI.VariableIndex, MOI.GreaterThan{T}),
        (MOI.VariableIndex, MOI.LessThan{T}),
    ]
end

function MOI.get(::SplitIntervalBridge, ::MOI.NumberOfVariables)::Int64
    return 2
end

function MOI.get(bridge::SplitIntervalBridge, ::MOI.ListOfVariableIndices)
    return [bridge.xp, bridge.xn]
end

# The constrained variables create VariableIndex-in-GreaterThan/LessThan
# constraints that must be reported as part of the bridge.

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{MOI.VariableIndex,MOI.GreaterThan{T}},
)::Int64 where {T}
    return 1 + (bridge.xn_lower !== nothing)
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{T}},
) where {T}
    ret = [
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{T}}(
            bridge.xp.value,
        ),
    ]
    if bridge.xn_lower !== nothing
        push!(ret, bridge.xn_lower)
    end
    return ret
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{MOI.VariableIndex,MOI.LessThan{T}},
)::Int64 where {T}
    return 1 + (bridge.xp_upper !== nothing)
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{T}},
) where {T}
    ret = [
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{T}}(bridge.xn.value),
    ]
    if bridge.xp_upper !== nothing
        push!(ret, bridge.xp_upper)
    end
    return ret
end

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
) where {T}
    return [bridge.lower]
end

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
) where {T}
    return [bridge.upper]
end

function MOI.get(
    ::SplitIntervalBridge{T,G},
    ::MOI.NumberOfConstraints{G,MOI.EqualTo{T}},
)::Int64 where {T,G}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T,G},
    ::MOI.ListOfConstraintIndices{G,MOI.EqualTo{T}},
) where {T,G}
    return [bridge.equality]
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::SplitIntervalBridge,
)
    return bridge.func
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::SplitIntervalBridge,
)
    return bridge.set
end

function MOI.delete(model::MOI.ModelLike, bridge::SplitIntervalBridge)
    MOI.delete(model, bridge.lower)
    MOI.delete(model, bridge.upper)
    MOI.delete(model, bridge.equality)
    MOI.delete(model, bridge.xp)
    MOI.delete(model, bridge.xn)
    return
end
