struct VerticalBridge <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements}
    equalities::Vector{MOI.ConstraintIndex}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{VerticalBridge},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::MOI.Complements,
)
    return VerticalBridge(reformulate_to_vertical!(model, func, set)...)
end

#=
    Parser for JuMP problems with complementarity constraints.
=#

_is_single_variable(::MOI.AbstractScalarFunction) = false
function _is_single_variable(func::MOI.ScalarAffineFunction)
    return length(func.terms) == 1 &&
           func.terms[1].coefficient == 1.0 &&
           iszero(func.constant)
end
function _is_single_variable(func::MOI.ScalarQuadraticFunction)
    if (
        length(func.quadratic_terms) == 0 &&
        length(func.affine_terms) == 1 &&
        func.affine_terms[1].coefficient == 1.0 &&
        iszero(func.constant)
    )
        return true
    else
        false
    end
end
function _is_single_variable(func::MOI.ScalarNonlinearFunction)
    return func.head == :+ && length(func.args) == 1 && isa(func.args[1], MOI.VariableIndex)
end
_get_variable(func::MOI.ScalarAffineFunction) = func.terms[1].variable
_get_variable(func::MOI.ScalarQuadraticFunction) = func.affine_terms[1].variable
_get_variable(func::MOI.ScalarNonlinearFunction) = func.args[1]


# TODO: add support for ScalarNonlinearTerm
function _parse_complementarity_constraint(fun::MOI.AbstractVectorFunction, n_comp)
    exprs = MOIU.scalarize(fun)
    @assert length(exprs) == 2*n_comp

    cc_lhs = MOI.AbstractScalarFunction[]
    cc_rhs = MOI.VariableIndex[]

    for i = 1:n_comp
        # Parse LHS
        t1 = exprs[i]
        t2 = exprs[i+n_comp]
        if _is_single_variable(t1)
            push!(cc_lhs, _get_variable(t1))
        else
            push!(cc_lhs, t1)
        end

        # Parse RHS
        isvar2 = _is_single_variable(t2)
        if !isvar2
            # The RHS should be a variable if we follow MOI's specs
            # TODO: we should decide if we should add support complementarity
            # between expressions (see Issue #2)
            error(
                "Right-hand-side should be a single variable in complementarity constraints.",
            )
        end
        push!(cc_rhs, _get_variable(t2))
    end

    return cc_lhs, cc_rhs
end

"""
    reformulate_to_vertical!(model::MOI.ModelLike)

Factorize all the complementarity constraints in `model` and formulate
an equivalent model in vertical form. The complementarity constraints involving
expressions are rewritten with a slack.

Once reformulated, the complementarity constraints involve only single variables.

Return a `Tuple{Vector{MOI.VariableIndex}, Vector{MOI.VariableIndex}}` storing the indexes
of the variables in the left-hand-side and the right-hand-side of each complementarity.

"""
function reformulate_to_vertical!(model::MOI.ModelLike, fun, set)
    equalities = MOI.ConstraintIndex[]
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]
    n_comp = div(set.dimension, 2)
    @assert !(fun isa MOI.VectorOfVariables)
    # Read each complementarity constraint and get corresponding indices
    cc_lhs, cc_rhs = _parse_complementarity_constraint(fun, n_comp)
    for (lhs, x2) in zip(cc_lhs, cc_rhs)
        # Check if x2 is bounded.
        lb, ub = MOIU.get_bounds(model, Float64, x2)
        if isinf(lb) && isinf(ub)
            # If x2 is unbounded, the LHS is directly converted to an equality constraint.
            push!(equalities, MOI.add_constraint(model, lhs, MOI.EqualTo{Float64}(0)))
        elseif isa(lhs, MOI.VariableIndex)
            # If lhs is a variable, no need to reformulate the
            # complementarity constraint using a slack.
            # TODO: we should check if the variable lhs is bounded.
            push!(ind_cc1, lhs)
            push!(ind_cc2, x2)
        else
            # Else, reformulate LHS using vertical form
            x1 = MOI.add_variable(model)
            new_lhs = MOIU.operate!(-, Float64, lhs, x1)
            push!(equalities, MOI.add_constraint(model, new_lhs, MOI.EqualTo{Float64}(0)))
            push!(ind_cc1, x1)
            push!(ind_cc2, x2)
        end
    end
    n_cc = length(ind_cc1)
    comp = MOI.VectorOfVariables([ind_cc1; ind_cc2])
    ci = MOI.add_constraint(model, comp, MOI.Complements(2*n_cc))
    return ci, equalities
end
