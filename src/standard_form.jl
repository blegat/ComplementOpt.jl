
"""
    reformulate_to_standard_form!(model::ModelLike)

The mixed-complementarity constraints ``x1 ⟂ (lb <= x2 <= ub)`` is rewriten with slack variables as
```
x2p = x2 - lb
x2n = ub - x2
x1 = x1p - x1n
0 <= x1p ⟂ x2p >= 0
0 <= x1n ⟂ x2n >= 0

```

The function assumes the model is passed in vertical form, and returns an error otherwise.
"""
function reformulate_to_standard_form!(model::MOI.ModelLike)
    if !is_vertical(model)
        error("Complementarity constraints should be reformulated in vertical form before reformulating them in standard form")
    end
    cc_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}())[1]
    fun = MOI.get(model, MOI.ConstraintFunction(), cc_cons)
    set = MOI.get(model, MOI.ConstraintSet(), cc_cons)
    MOI.delete(model, cc_cons)
    n_comp = div(set.dimension, 2)

    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]

    for cc in 1:n_comp
        x1 = fun.variables[cc]
        x2 = fun.variables[cc + n_comp]
        # Get bounds on x1
        lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
        # Get bounds on x2
        lb2, ub2 = MOIU.get_bounds(model, Float64, x2)

        # If lb2 = ub2, the left-hand-side x1 is free and we discard this
        # complementarity constraint.
        if lb2 == ub2
            continue
        end

        if isfinite(lb2) && isinf(ub2)
            # Ensure x1 is well posed
            if (isfinite(lb1) && !iszero(lb1)) || isfinite(ub1)
                error("The problem does not follow MOI convention for mixed-complementarity constraints: the LHS lower bound is not zero.")
            end
            if isinf(lb1)
                MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
            end
            # Add a slack x2p = x2 - lb2 if lb2 ≠ 0
            if !iszero(lb2)
                x2p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    1.0 * x2 - 1.0 * x2p,
                    MOI.EqualTo(lb2),
                )
                push!(ind_cc1, x1)
                push!(ind_cc2, x2p)
            else
                # Nothing to change
                push!(ind_cc1, x1)
                push!(ind_cc2, x2)
            end
        elseif isinf(lb2) && isfinite(ub2)
            # Ensure x1 is well posed
            if (isfinite(ub1) && !iszero(ub1)) || isfinite(lb1)
                error("The problem does not follow MOI convention for mixed-complementarity constraints: the RHS lower bound is not zero.")
            end
            # Add a slack x1n = -x1
            x1n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            # x1n = -x1
            MOI.add_constraint(
                model,
                1.0 * x1 + 1.0 * x1n,
                MOI.EqualTo(0.0),
            )
            # Add a slack x2n = ub2 - x2
            x2n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                1.0*x2 + 1.0*x2n,
                MOI.EqualTo(ub2),
            )
            push!(ind_cc1, x1n)
            push!(ind_cc2, x2n)
        elseif isfinite(lb2) && isfinite(ub2)
            # Reformulate the mixed-complementarity constraint as two complementarity constraints.
            # First, ensure x1 is well posed
            if isfinite(lb1) || isfinite(ub1)
                error("The problem does not follow MOI convention for mixed-complementarity constraints: the LHS should be free.")
            end

            # Add two slacks such that x1 = x1p - x1n
            x1p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x1n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                1.0*x1 - 1.0*x1p + 1.0*x1n,
                MOI.EqualTo(0.0),
            )
            push!(ind_cc1, x1p)
            if !iszero(lb2)
                x2p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    1.0*x2 - 1.0*x2p,
                    MOI.EqualTo(lb2),
                )
                push!(ind_cc2, x2p)
            else
                push!(ind_cc2, x2)
            end

            x2n, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                1.0*x2 + 1.0*x2n,
                MOI.EqualTo(ub2),
            )
            push!(ind_cc1, x1n)
            push!(ind_cc2, x2n)
        else
            error("Problem is not well specified")
        end
    end

    @assert length(ind_cc1) == length(ind_cc2)
    n_cc = length(ind_cc1)
    comp = MOI.VectorOfVariables([ind_cc1; ind_cc2])
    MOI.add_constraint(model, comp, MOI.Complements(2*n_cc))

    return ind_cc1, ind_cc2
end

