struct NonlinearReformulation end

"""
    reformulate_as_nonlinear_program!(model::MOI.ModelLike; relaxation=0.0)

Reformulate a model with complementarity constraints written in vertical form
to an equivalent (degenerate) nonlinear program.

If the complementarity constraints are not in vertical form, an error is thrown.

"""
function reformulate_as_nonlinear_program!(model::MOI.ModelLike; relaxation=0.0)
    if !is_vertical(model)
        error("Complementarity constraints should be reformulated in vertical form before applying nonlinear reformulation")
    end

    cc_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}())[1]
    fun = MOI.get(model, MOI.ConstraintFunction(), cc_cons)
    set = MOI.get(model, MOI.ConstraintSet(), cc_cons)
    n_comp = div(set.dimension, 2)

    ind_cc = []
    for cc in 1:n_comp
        x1 = fun.variables[cc]
        x2 = fun.variables[cc + n_comp]
        # Get bounds on x1
        lb1, ub1 = MOIU.get_bounds(model, Float64, x2)
        # Get bounds on x2
        lb2, ub2 = MOIU.get_bounds(model, Float64, x2)
        # x2 should have at least one bound, otherwise x1 becomes a fixed variable
        @assert isinf(lb2) || isinf(ub2)
        quad_terms = [
            MOI.ScalarQuadraticTerm(1.0, x1, x2)
        ]
        if !isinf(lb2) && isinf(ub2)
            # x1 ⟂ (lb <= x2)   ≡  0 <= x1 ; lb <= x2 ; x1 ( x2 - lb) <= 0
            if isinf(lb1)
                MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
            else
                @assert lb1 == 0.0 # ensure we follow MOI's convention
                # TODO: what should we do if ub1 is finite?
            end
            idc = MOI.add_constraint(
                model,
                MOI.ScalarQuadraticFunction(
                    quad_terms,
                    [MOI.ScalarAffineTerm(-lb2, x1)],
                    0.0,
                ),
                MOI.LessThan(relaxation),
            )
            push!(ind_cc, idc)
        elseif isinf(lb2) && !isinf(ub2)
            # x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; x1 ( x2 - ub) <= 0
            if isinf(ub1)
                MOI.add_constraint(model, x1, MOI.LessThan(0.0))
            else
                @assert ub1 == 0.0 # ensure we follow MOI's convention
                # TODO: what should we do if lb1 is finite?
            end
            idc = MOI.add_constraint(
                model,
                MOI.ScalarQuadraticFunction(
                    quad_terms,
                    [MOI.ScalarAffineTerm(-ub2, x1)],
                    0.0,
                ),
                MOI.LessThan(relaxation),
            )
            push!(ind_cc, idc)
        else
            # x1 ⟂ (lb <= x2 <= ub) ≡  lb <= x2 <= ub ; x1 (x2 - lb)  <= 0 ; x1 (x2 - ub) <= 0
            idc = MOI.add_constraint(
                model,
                MOI.ScalarQuadraticFunction(
                    quad_terms,
                    [MOI.ScalarAffineTerm(-lb2, x1)],
                    0.0,
                ),
                MOI.LessThan(relaxation),
            )
            push!(ind_cc, idc)
            idc = MOI.add_constraint(
                model,
                MOI.ScalarQuadraticFunction(
                    quad_terms,
                    [MOI.ScalarAffineTerm(-ub2, x1)],
                    0.0,
                ),
                MOI.LessThan(relaxation),
            )
            push!(ind_cc, idc)
        end
    end
    MOI.delete(model, cc_cons)
    return ind_cc
end
