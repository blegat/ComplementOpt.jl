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
[`MOI.Bridges.final_touch`](@ref), if the activity `x` has finite lower
and upper bounds `[lb, ub]`, the bounds `xp <= max(ub, 0)` and
`xn >= min(lb, 0)` are added. This is done in
[`MOI.Bridges.final_touch`](@ref) and not in
`MOI.Bridges.Constraint.bridge_constraint` because the bounds of `x` may
be set after the constraint is bridged. These bounds are needed by
[`MOI.Bridges.Constraint.SOS1ToMILPBridge`](@ref) in case the inner solver
does not support [`MOI.SOS1`](@ref) constraints. Since that bridge
requires both bounds anyway (`xp` and `xn` each appear in an
[`MOI.SOS1`](@ref) constraint), we use the same
`MOI.Utilities.get_bounds` method and add nothing when either bound is
infinite.

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
    # `lower` and `upper` are created in `final_touch` (after `xp`/`xn` are
    # bounded) so that the bridges they trigger (e.g. `VerticalBridge`) see the
    # finite bounds. They are `nothing` until then.
    lower::Union{
        Nothing,
        MOI.ConstraintIndex{
            MOI.VectorOfVariables,
            ComplementsWithSetType{MOI.GreaterThan{T}},
        },
    }
    upper::Union{
        Nothing,
        MOI.ConstraintIndex{
            MOI.VectorOfVariables,
            ComplementsWithSetType{MOI.LessThan{T}},
        },
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
    reformulation::Union{Nothing,AbstractComplementarityRelaxation}
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
        # Extract the variable from a function wrapping a single variable
        @assert _is_single_variable(y)
        _get_variable(y)
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
    # `lower` and `upper` are created in `final_touch`, once `xp`/`xn` are bounded.
    return SplitIntervalBridge{T,G,typeof(func)}(
        nothing,
        nothing,
        equality,
        xp,
        xn,
        func,
        set,
        nothing,
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
    bridge.reformulation = value
    if bridge.lower !== nothing
        MOI.set(model, attr, bridge.lower, value)
    end
    if bridge.upper !== nothing
        MOI.set(model, attr, bridge.upper, value)
    end
    return
end

# The bounds of the activity `x`, or `nothing` if they cannot be computed or the
# domain is not bounded. `get_bounds` is only defined for affine functions and
# variables.
function _activity_bounds(model, ::Type{T}, x) where {T}
    if x isa MOI.VariableIndex || x isa MOI.ScalarAffineFunction{T}
        cache = Dict{MOI.VariableIndex,NTuple{2,T}}()
        return MOI.Utilities.get_bounds(model, cache, x)
    end
    return nothing
end

MOI.Bridges.needs_final_touch(::SplitIntervalBridge) = true

function MOI.Bridges.final_touch(
    bridge::SplitIntervalBridge{T},
    model::MOI.ModelLike,
) where {T}
    # Bound `xp` and `xn` from the activity bounds. This is done before creating
    # `lower`/`upper` below so that the bridges they trigger (e.g.
    # `VerticalBridge`) can compute the bounds of the slacks they introduce.
    ret = _activity_bounds(model, T, first(MOI.Utilities.eachscalar(bridge.func)))
    if ret === nothing
        if bridge.xp_upper !== nothing
            MOI.delete(model, bridge.xp_upper)
            bridge.xp_upper = nothing
        end
        if bridge.xn_lower !== nothing
            MOI.delete(model, bridge.xn_lower)
            bridge.xn_lower = nothing
        end
    else
        lb, ub = ret
        upper = MOI.LessThan(max(ub, zero(T)))
        if bridge.xp_upper === nothing
            bridge.xp_upper = MOI.add_constraint(model, bridge.xp, upper)
        elseif MOI.get(model, MOI.ConstraintSet(), bridge.xp_upper) != upper
            MOI.set(model, MOI.ConstraintSet(), bridge.xp_upper, upper)
        end
        lower = MOI.GreaterThan(min(lb, zero(T)))
        if bridge.xn_lower === nothing
            bridge.xn_lower = MOI.add_constraint(model, bridge.xn, lower)
        elseif MOI.get(model, MOI.ConstraintSet(), bridge.xn_lower) != lower
            MOI.set(model, MOI.ConstraintSet(), bridge.xn_lower, lower)
        end
    end
    # Create the two complementarity constraints (only once).
    if bridge.lower === nothing
        y_var = MOI.Utilities.scalarize(bridge.func)[2]
        if !(y_var isa MOI.VariableIndex)
            y_var = _get_variable(y_var)
        end
        bridge.lower = MOI.add_constraint(
            model,
            MOI.VectorOfVariables([bridge.xp, y_var]),
            ComplementsWithSetType{MOI.GreaterThan{T}}(2),
        )
        bridge.upper = MOI.add_constraint(
            model,
            MOI.VectorOfVariables([bridge.xn, y_var]),
            ComplementsWithSetType{MOI.LessThan{T}}(2),
        )
        if bridge.reformulation !== nothing
            MOI.set(
                model,
                ComplementarityReformulation(),
                bridge.lower,
                bridge.reformulation,
            )
            MOI.set(
                model,
                ComplementarityReformulation(),
                bridge.upper,
                bridge.reformulation,
            )
        end
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
    bridge::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
)::Int64 where {T}
    return bridge.lower === nothing ? 0 : 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
) where {T}
    F, S = MOI.VectorOfVariables, ComplementsWithSetType{MOI.GreaterThan{T}}
    return bridge.lower === nothing ? MOI.ConstraintIndex{F,S}[] :
           [bridge.lower]
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
)::Int64 where {T}
    return bridge.upper === nothing ? 0 : 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
) where {T}
    F, S = MOI.VectorOfVariables, ComplementsWithSetType{MOI.LessThan{T}}
    return bridge.upper === nothing ? MOI.ConstraintIndex{F,S}[] :
           [bridge.upper]
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
    if bridge.lower !== nothing
        MOI.delete(model, bridge.lower)
    end
    if bridge.upper !== nothing
        MOI.delete(model, bridge.upper)
    end
    MOI.delete(model, bridge.equality)
    # The bounds `xp_upper`/`xn_lower` are deleted together with `xp`/`xn`.
    MOI.delete(model, bridge.xp)
    MOI.delete(model, bridge.xn)
    return
end
