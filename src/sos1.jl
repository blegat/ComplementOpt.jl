"""
    SOS1Relaxation <: AbstractComplementarityRelaxation

Reformulate complementarity constraints using SOS1 constraints.
The complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
SOS1(a, b)
```
where the SOS1 constraint enforces that at most one of `a` and `b` is nonzero.

When the variables are shifted (e.g., `b ≥ lb` with `lb ≠ 0`), slack variables
are introduced so that the SOS1 constraint operates on nonneg variables that are
zero at the bound.

!!! note
    If the solver does not support SOS1 constraints natively, MOI will bridge
    them to MILP using the `SOS1ToMILPBridge`, which requires all variables in
    the SOS1 constraint to have finite bounds. In that case, make sure the inner
    solver is wrapped with `MOI.Bridges.full_bridge_optimizer` and all problem
    variables have finite bounds.

"""
struct SOS1Relaxation <: AbstractComplementarityRelaxation end

# x1 ⟂ (lb2 <= x2)   ≡  0 <= x1 ; lb2 <= x2 ; SOS1(x1, x2 - lb2)
function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    ::SOS1Relaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, _ = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0
    end
    # Create slack s = x2 - lb2 >= 0
    s, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    eq = MOI.add_constraint(model, 1.0 * x2 - 1.0 * s, MOI.EqualTo(lb2))
    sos = MOI.add_constraint(model, MOI.VectorOfVariables([x1, s]), MOI.SOS1([1.0, 2.0]))
    return [eq, sos]
end

# x1 ⟂ (x2 <= ub2)   ≡  x1 <= 0 ; x2 <= ub2 ; SOS1(-x1, ub2 - x2)
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    ::SOS1Relaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    _, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.LessThan(0.0))
    else
        @assert ub1 == 0.0
    end
    # s1 = -x1 >= 0
    s1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    eq1 = MOI.add_constraint(model, 1.0 * x1 + 1.0 * s1, MOI.EqualTo(0.0))
    # s2 = ub2 - x2 >= 0
    s2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    eq2 = MOI.add_constraint(model, 1.0 * x2 + 1.0 * s2, MOI.EqualTo(ub2))
    sos = MOI.add_constraint(model, MOI.VectorOfVariables([s1, s2]), MOI.SOS1([1.0, 2.0]))
    return [eq1, eq2, sos]
end

# x1 ⟂ (lb2 <= x2 <= ub2)
# Split x1 = p - q with p, q >= 0, then:
#   SOS1(p, x2 - lb2)  and  SOS1(q, ub2 - x2)
function _relax_complementarity_range!(
    model::MOI.ModelLike,
    ::SOS1Relaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    range2 = ub2 - lb2
    # Compute upper bounds for p and q from x1 bounds
    ub_p = isinf(ub1) ? Inf : max(ub1, 0.0)
    ub_q = isinf(lb1) ? Inf : max(-lb1, 0.0)
    # Split x1 = p - q, p >= 0, q >= 0
    p, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    if isfinite(ub_p)
        MOI.add_constraint(model, p, MOI.LessThan(ub_p))
    end
    q, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    if isfinite(ub_q)
        MOI.add_constraint(model, q, MOI.LessThan(ub_q))
    end
    eq_split = MOI.add_constraint(model, 1.0 * x1 - 1.0 * p + 1.0 * q, MOI.EqualTo(0.0))
    # s_lb = x2 - lb2 ∈ [0, range2]
    s_lb, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, s_lb, MOI.LessThan(range2))
    eq_lb = MOI.add_constraint(model, 1.0 * x2 - 1.0 * s_lb, MOI.EqualTo(lb2))
    # s_ub = ub2 - x2 ∈ [0, range2]
    s_ub, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, s_ub, MOI.LessThan(range2))
    eq_ub = MOI.add_constraint(model, 1.0 * x2 + 1.0 * s_ub, MOI.EqualTo(ub2))
    # SOS1 constraints
    sos1 = MOI.add_constraint(model, MOI.VectorOfVariables([p, s_lb]), MOI.SOS1([1.0, 2.0]))
    sos2 = MOI.add_constraint(model, MOI.VectorOfVariables([q, s_ub]), MOI.SOS1([1.0, 2.0]))
    return [eq_split, eq_lb, eq_ub, sos1, sos2]
end
