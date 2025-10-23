
const VF = Union{
    MOI.VectorAffineFunction,
    MOI.VectorQuadraticFunction,
}

function has_complementarity(model::MOI.ModelLike)
    return any(MOI.get(model, MOI.ListOfConstraintTypesPresent())) do (F, S)
        S == MOI.Complements
    end
end

#=
    Parser for JuMP problems with complementarity constraints.
=#

_is_single_variable(::MOI.AbstractScalarFunction) = false
function _is_single_variable(func::MOI.ScalarAffineFunction)
    return length(func.terms) == 1 && func.terms[1].coefficient == 1.0 && iszero(func.constant)
end
function _is_single_variable(func::MOI.ScalarQuadraticFunction)
    if (length(func.quadratic_terms) == 0
        && length(func.affine_terms) == 1
        && func.affine_terms[1].coefficient == 1.0
        && iszero(func.constant)
    )
        return true
    else
        false
    end
end
_get_variable(func::MOI.ScalarAffineFunction) = func.terms[1].variable
_get_variable(func::MOI.ScalarQuadraticFunction) = func.affine_terms[1].variable

# Reformulate
function _vertical_formulation!(model, terms::Vector{MOI.ScalarAffineTerm{T}}, c::T) where T
    x = MOI.add_variable(model)
    push!(terms, MOI.ScalarAffineTerm{T}(-one(T), x))
    func = MOI.ScalarAffineFunction{T}(terms, zero(T))
    MOI.add_constraint(model, func, MOI.EqualTo{T}(-c))
    return x
end

function _add_slack!(model, func::MOI.ScalarAffineFunction{T}) where T
    x = MOI.add_variable(model)
    push!(func.terms, MOI.ScalarAffineTerm{T}(-one(T), x))
    return x
end
function _add_slack!(model, func::MOI.ScalarQuadraticFunction{T}) where T
    x = MOI.add_variable(model)
    push!(func.affine_terms, MOI.ScalarAffineTerm{T}(-one(T), x))
    return x
end

# TODO: add support for ScalarNonlinearTerm
function _parse_complementarity_constraint(fun::MOI.AbstractVectorFunction, n_comp)
    exprs = MOIU.scalarize(fun)
    @assert length(exprs) == 2*n_comp

    cc_lhs = MOI.AbstractScalarFunction[]
    cc_rhs = MOI.VariableIndex[]

    for i in 1:n_comp
        # Parse LHS
        t1 = exprs[i]
        t2 = exprs[i + n_comp]
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
            error("Right-hand-side should be a single variable in complementarity constraints.")
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
function reformulate_to_vertical!(model::MOI.ModelLike)
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]

    contypes = MOI.get(model, MOI.ListOfConstraintTypesPresent())
    for (F, S) in contypes
        # Parse only complementarity constraints
        if S == MOI.Complements
            conindices = MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
            for cidx in conindices
                fun = MOI.get(model, MOI.ConstraintFunction(), cidx)
                set = MOI.get(model, MOI.ConstraintSet(), cidx)
                n_comp = div(set.dimension, 2)
                if isa(fun, MOI.VectorOfVariables)
                    append!(ind_cc1, fun.variables[1:n_comp])
                    append!(ind_cc2, fun.variables[n_comp+1:end])
                elseif isa(fun, VF)
                    # Read each complementarity constraint and get corresponding indices
                    cc_lhs, cc_rhs = _parse_complementarity_constraint(fun, n_comp)
                    for (lhs, x2) in zip(cc_lhs, cc_rhs)
                        # Check if x2 is bounded.
                        lb, ub = MOIU.get_bounds(model, Float64, x2)
                        if isinf(lb) && isinf(ub)
                            # If x2 is unbounded, the LHS is directly converted to an equality constraint.
                            MOI.add_constraint(model, lhs, MOI.EqualTo{Float64}(0))
                        elseif isa(lhs, MOI.VariableIndex)
                            # If lhs is a variable, no need to reformulate the
                            # complementarity constraint using a slack.
                            # TODO: we should check if the variable lhs is bounded.
                            push!(ind_cc1, lhs)
                            push!(ind_cc2, x2)
                        else
                            # Else, reformulate LHS using vertical form
                            x1 = _add_slack!(model, lhs)
                            MOI.add_constraint(model, lhs, MOI.EqualTo{Float64}(0))
                            push!(ind_cc1, x1)
                            push!(ind_cc2, x2)
                        end
                    end
                else
                    error("Complementary constraints formulated with $(typeof(fun)) are not yet supported")
                end
                # We delete the complementarity constraints
                MOI.delete(model, cidx)
            end
        end
    end
    n_cc = length(ind_cc1)
    comp = MOI.VectorOfVariables([ind_cc1; ind_cc2])
    MOI.add_constraint(model, comp, MOI.Complements(2*n_cc))
    return ind_cc1, ind_cc2
end

function is_vertical(model::MOI.ModelLike)
    one_vov = false
    for (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        if S == MOI.Complements
            if F == MOI.VectorOfVariables
                if MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 1
                    one_vov = true
                else
                    return false
                end
            else
                return false
            end
        end
    end
    return one_vov
end
