
"""
    AbstractComplementarityRelaxation

Abstract type to implement any complementarity function ``\\psi``.

"""
abstract type AbstractComplementarityRelaxation end

"""
    reformulate_as_nonlinear_program!(model::MOI.ModelLike, relaxation::AbstractComplementarityRelaxation)

Reformulate a MOI model with complementarity constraints written in vertical form
as a nonlinear program. The second argument specifies the relaxation used to reformulate
the complementarity constraints.

If the complementarity constraints are not in vertical form, an error is thrown.

"""
function reformulate_as_nonlinear_program!(
    model::MOI.ModelLike,
    relaxation::AbstractComplementarityRelaxation,
)
    if !is_vertical(model)
        error(
            "Complementarity constraints should be reformulated in vertical form before applying nonlinear reformulation",
        )
    end

    cc_cons = MOI.get(
        model,
        MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.Complements}(),
    )[1]
    fun = MOI.get(model, MOI.ConstraintFunction(), cc_cons)
    set = MOI.get(model, MOI.ConstraintSet(), cc_cons)
    n_comp = div(set.dimension, 2)

    ind_cc = []
    for cc = 1:n_comp
        x1 = fun.variables[cc]
        x2 = fun.variables[cc+n_comp]
        # Get bounds on x2
        lb2, ub2 = MOIU.get_bounds(model, Float64, x2)
        # x2 should have at least one bound, otherwise x1 becomes a fixed variable
        @assert isfinite(lb2) || isfinite(ub2)
        if !isinf(lb2) && isinf(ub2)
            idc = _relax_complementarity_lower_bound!(model, relaxation, x1, x2, lb2, ub2)
        elseif isinf(lb2) && !isinf(ub2)
            idc = _relax_complementarity_upper_bound!(model, relaxation, x1, x2, lb2, ub2)
        else
            idc = _relax_complementarity_range!(model, relaxation, x1, x2, lb2, ub2)
        end
        append!(ind_cc, idc)
    end
    MOI.delete(model, cc_cons)
    return ind_cc
end


"""
    ScholtesRelaxation <: AbstractComplementarityRelaxation

Implement the Scholtes relaxation.
For `tau ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
a . b ≤ tau
```

If `tau` is equal to 0, then the relaxation is exact.

"""
struct ScholtesRelaxation{T} <: AbstractComplementarityRelaxation
    tau::T
end

# x1 ⟂ (lb <= x2)   ≡  0 <= x1 ; lb <= x2 ; x1 ( x2 - lb) <= 0
function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.tau))
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; x1 ( x2 - ub) <= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.tau))
    return [idc]
end

# x1 ⟂ (lb <= x2 <= ub) ≡  lb <= x2 <= ub ; x1 (x2 - lb)  <= 0 ; x1 (x2 - ub) <= 0
function _relax_complementarity_range!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    idc1 = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.tau))
    idc2 = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.tau))
    return [idc1, idc2]
end

"""
    FischerBurmeisterRelaxation <: AbstractComplementarityRelaxation

Implement the Fischer-Burmeister relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
a + b - sqrt((a + b)^2 + epsilon) ≤ 0
```

"""
struct FischerBurmeisterRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _min_eps(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :-,
        Any[
            1.0*a+1.0*b,
            MOI.ScalarNonlinearFunction(:sqrt, Any[(1.0*a)^2+(1.0*b)^2+eps^2]),
        ],
    )
end

function _max_eps(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :+,
        Any[
            1.0*a+1.0*b,
            MOI.ScalarNonlinearFunction(:sqrt, Any[(1.0*a)^2+(1.0*b)^2+eps^2]),
        ],
    )
end

# x1 ⟂ (lb <= x2)   ≡  0 <= x1 ; lb <= x2 ; min(x1, x2 - lb) <= 0
function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _min_eps(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _max_eps(x1, x2 - ub2, relaxation.epsilon),
        MOI.GreaterThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (lb <= x2 <= ub) ≡  lb <= x2 <= ub ; min(x1, x2 - lb) <= 0 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_range!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    idc1 = MOI.add_constraint(
        model,
        _min_eps(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    idc2 = MOI.add_constraint(
        model,
        _max_eps(x1, x2 - ub2, relaxation.epsilon),
        MOI.GreaterThan(0.0),
    )
    return [idc1, idc2]
end


"""
    LiuFukushimaRelaxation <: AbstractComplementarityRelaxation

Implement the Liu-Fukushima relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
a . b ≤ epsilon^2
(a + epsilon) . (b + epsilon) ≥ epsilon^2

```
"""
struct LiuFukushimaRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::LiuFukushimaRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    idc1 = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.epsilon^2))
    idc2 = MOI.add_constraint(
        model,
        (x1 + relaxation.epsilon) * (x2 - lb2 + relaxation.epsilon),
        MOI.GreaterThan(relaxation.epsilon^2),
    )

    # Remove bounds
    _remove_bounds!(model, x1)
    cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2.value)
    MOI.delete(model, cidx)

    return [idc1, idc2]
end

function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::LiuFukushimaRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    idc1 = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.epsilon^2))
    idc2 = MOI.add_constraint(
        model,
        (x1 - relaxation.epsilon) * (x2 - ub2 - relaxation.epsilon),
        MOI.GreaterThan(relaxation.epsilon^2),
    )

    # Remove bounds
    _remove_bounds!(model, x1)
    cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2.value)
    MOI.delete(model, cidx)

    return [idc1, idc2]
end


"""
    KanzowSchwarzRelaxation <: AbstractComplementarityRelaxation

Implement the Kanzow-Schwarz relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
ϕ(a, b) ≤ 0

```
with the function `ϕ`:
```
ϕ(a, b) = (a - epsilon) . (b - epsilon)                 if a + b > 2 epsilon
          -0.5 ((a -epsilon)^2 + (b - epsilon)^2)       otherwise

```

"""
struct KanzowSchwarzRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _kanzow_schwarz_relaxation(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :ifelse,
        Any[
            MOI.ScalarNonlinearFunction(:>, Any[1.0*a+1.0*b, 2*eps]),
            (a-eps)*(b-eps),
            - 0.5*((a-eps)^2+(b-eps)^2),
        ],
    )
end

function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::KanzowSchwarzRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _kanzow_schwarz_relaxation(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::KanzowSchwarzRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _kanzow_schwarz_relaxation(-x1, ub2 - x2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end
