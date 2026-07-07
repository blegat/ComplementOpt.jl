# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    VerticalBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge

`VerticalBridge` implements the following reformulation:

  * `f(x)` in `S` into `[x₁, x₂]` in `S`

where expression-based complementarity constraints are converted to vertical
form by introducing slack variables. If the left-hand side is an expression,
a slack variable `x₁` is created with an equality `lhs = x₁`. If the
right-hand side variable `x₂` is unbounded, the left-hand side is converted
to an equality constraint instead.

## Source node

`VerticalBridge` supports:

  * [`MOI.AbstractVectorFunction`](@ref) in [`MOI.Complements`](@ref)
  * [`MOI.AbstractVectorFunction`](@ref) in [`ComplementsWithSetType{S}`](@ref)

## Target nodes

`VerticalBridge` creates:

  * [`MOI.VectorOfVariables`](@ref) in `S`
  * [`MOI.ScalarAffineFunction{T}`](@ref) in [`MOI.EqualTo{T}`](@ref)

"""
struct VerticalBridge{T,S<:MOI.AbstractVectorSet} <:
       MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    equalities::Vector{MOI.ConstraintIndex}
    slacks::Vector{MOI.VariableIndex}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{VerticalBridge{T,MOI.Complements}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::MOI.Complements,
) where {T}
    ci, equalities, slacks = reformulate_to_vertical!(model, T, func, set)
    return VerticalBridge{T,MOI.Complements}(ci, equalities, slacks)
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{VerticalBridge{T,ComplementsWithSetType{S}}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::ComplementsWithSetType{S},
) where {T,S}
    ci, equalities, slacks = reformulate_to_vertical!(model, T, func, set)
    return VerticalBridge{T,ComplementsWithSetType{S}}(ci, equalities, slacks)
end

function MOI.supports_constraint(
    ::Type{<:VerticalBridge},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
)
    return true
end

function MOI.supports_constraint(
    ::Type{<:VerticalBridge},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:VerticalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) where {T}
    return VerticalBridge{T,MOI.Complements}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:VerticalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    return VerticalBridge{T,ComplementsWithSetType{S}}
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:VerticalBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{VerticalBridge{T,S}},
) where {T,S}
    return Tuple{Type,Type}[
        (MOI.VectorOfVariables, S),
        (MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}),
    ]
end

function MOI.get(bridge::VerticalBridge, ::MOI.NumberOfVariables)::Int64
    return length(bridge.slacks)
end

function MOI.get(bridge::VerticalBridge, ::MOI.ListOfVariableIndices)
    return copy(bridge.slacks)
end

function MOI.get(
    ::VerticalBridge{T,S},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
)::Int64 where {T,S}
    return 1
end

function MOI.get(
    bridge::VerticalBridge{T,S},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {T,S}
    return [bridge.constraint]
end

function MOI.get(
    bridge::VerticalBridge{T},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
)::Int64 where {T}
    return count(
        ci ->
            ci isa
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
        bridge.equalities,
    )
end

function MOI.get(
    bridge::VerticalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
) where {T}
    return MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}}[
        ci for ci in bridge.equalities if
        ci isa MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}}
    ]
end

function MOI.delete(model::MOI.ModelLike, bridge::VerticalBridge)
    MOI.delete(model, bridge.constraint)
    for ci in bridge.equalities
        MOI.delete(model, ci)
    end
    for vi in bridge.slacks
        MOI.delete(model, vi)
    end
    return
end

function MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:VerticalBridge},
)
    return true
end

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::VerticalBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.constraint, value)
    return
end

#=
    Parser for JuMP problems with complementarity constraints.
=#

function _is_single_variable(func::MOI.ScalarAffineFunction)
    return length(func.terms) == 1 &&
           func.terms[1].coefficient == 1.0 &&
           iszero(func.constant)
end

function _is_single_variable(func::MOI.ScalarQuadraticFunction)
    return (
        length(func.quadratic_terms) == 0 &&
        length(func.affine_terms) == 1 &&
        func.affine_terms[1].coefficient == 1.0 &&
        iszero(func.constant)
    )
end

function _is_single_variable(func::MOI.ScalarNonlinearFunction)
    return func.head == :+ &&
           length(func.args) == 1 &&
           isa(func.args[1], MOI.VariableIndex)
end

_get_variable(func::MOI.ScalarAffineFunction) = func.terms[1].variable

_get_variable(func::MOI.ScalarQuadraticFunction) = func.affine_terms[1].variable

_get_variable(func::MOI.ScalarNonlinearFunction) = func.args[1]

# A single variable shifted by a constant, e.g. `x - 1`. Such a right-hand side
# is produced by `ComplementsVectorizeBridge` when it shifts the second variable
# by a nonzero bound.
function _is_shifted_variable(func::MOI.ScalarAffineFunction)
    return length(func.terms) == 1 && isone(func.terms[1].coefficient)
end

_is_shifted_variable(::MOI.AbstractScalarFunction) = false

# TODO: add support for ScalarNonlinearTerm
function _parse_complementarity_constraint(
    fun::MOI.AbstractVectorFunction,
    n_comp,
)
    exprs = MOI.Utilities.scalarize(fun)
    @assert length(exprs) == 2*n_comp
    cc_lhs = MOI.AbstractScalarFunction[exprs[i] for i in 1:n_comp]
    cc_rhs = MOI.AbstractScalarFunction[exprs[i+n_comp] for i in 1:n_comp]
    return cc_lhs, cc_rhs
end

# Return a single variable representing `expr`. If `expr` is already a single
# variable, it is returned directly. Otherwise, a slack variable `s` is created
# together with the equality `expr == s` and `s` is returned.
function _to_variable!(model, ::Type{T}, expr, slacks, equalities) where {T}
    if _is_single_variable(expr)
        return _get_variable(expr)
    end
    s = MOI.add_variable(model)
    push!(slacks, s)
    new_expr = MOI.Utilities.operate(-, T, expr, s)
    push!(
        equalities,
        MOI.Utilities.normalize_and_add_constraint(
            model,
            new_expr,
            MOI.EqualTo{T}(zero(T)),
        ),
    )
    return s
end

# Resolve the right-hand side of a complementarity constraint to a single
# variable. Following MOI's specs, it should be a single variable, but it may be
# a single variable shifted by a constant (introduced by
# `ComplementsVectorizeBridge`), in which case a slack is added. Genuine
# multi-variable expressions are not supported (see Issue #2).
function _rhs_variable!(model, ::Type{T}, expr, slacks, equalities) where {T}
    if !_is_single_variable(expr) && !_is_shifted_variable(expr)
        error(
            "Right-hand-side should be a single variable in complementarity constraints.",
        )
    end
    return _to_variable!(model, T, expr, slacks, equalities)
end

"""
    reformulate_to_vertical!(model::MOI.ModelLike, ::Type{T}, fun, set)

Factorize all the complementarity constraints in `model` and formulate
an equivalent model in vertical form. The complementarity constraints involving
expressions are rewritten with a slack. `T` is the coefficient type used for
the generated equality constraints.

Once reformulated, the complementarity constraints involve only single variables.
"""
function reformulate_to_vertical!(
    model::MOI.ModelLike,
    ::Type{T},
    fun,
    set,
) where {T}
    equalities = MOI.ConstraintIndex[]
    slacks = MOI.VariableIndex[]
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]
    n_comp = div(set.dimension, 2)
    @assert !(fun isa MOI.VectorOfVariables)
    # Read each complementarity constraint and get corresponding indices
    cc_lhs, cc_rhs = _parse_complementarity_constraint(fun, n_comp)
    for (lhs, rhs) in zip(cc_lhs, cc_rhs)
        x2 = _rhs_variable!(model, T, rhs, slacks, equalities)
        if set isa MOI.Complements
            # Check if x2 is bounded.
            lb, ub = MOI.Utilities.get_bounds(model, T, x2)
            if isinf(lb) && isinf(ub)
                # If x2 is unbounded, the LHS is directly converted to an equality constraint.
                push!(
                    equalities,
                    MOI.Utilities.normalize_and_add_constraint(
                        model,
                        lhs,
                        MOI.EqualTo{T}(zero(T)),
                    ),
                )
                continue
            end
        end
        # If lhs is a variable, no need to reformulate the complementarity
        # constraint using a slack.
        # TODO: we should check if the variable lhs is bounded.
        x1 = _to_variable!(model, T, lhs, slacks, equalities)
        push!(ind_cc1, x1)
        push!(ind_cc2, x2)
    end
    n_cc = length(ind_cc1)
    comp = MOI.VectorOfVariables([ind_cc1; ind_cc2])
    S = typeof(set)
    if set isa MOI.Complements
        ci = MOI.add_constraint(model, comp, MOI.Complements(2*n_cc))
    else
        ci = MOI.add_constraint(model, comp, S(2*n_cc))
    end
    return ci, equalities, slacks
end
