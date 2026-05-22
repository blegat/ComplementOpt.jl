# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    NonlinearBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge

`NonlinearBridge` implements the following reformulation:

  * `[x₁, x₂]` in [`ComplementsWithSetType{S}`](@ref) into nonlinear
    constraints using a complementarity relaxation method

The relaxation method is determined by the
[`AbstractComplementarityRelaxation`](@ref) attribute (default:
[`ScholtesRelaxation`](@ref)). The set type `S` determines the bound case:

  * `S <: Union{MOI.Nonnegatives, MOI.GreaterThan}`: lower-bound relaxation
  * `S <: Union{MOI.Nonpositives, MOI.LessThan}`: upper-bound relaxation
  * `S <: MOI.Interval`: range relaxation
  * `S == MOI.Zeros`: equality (trivial complementarity)

## Source node

`NonlinearBridge` supports:

  * [`MOI.VectorOfVariables`](@ref) in [`ComplementsWithSetType{S}`](@ref)

## Target nodes

`NonlinearBridge` creates nonlinear constraints depending on the relaxation
method (for example, quadratic or `ScalarNonlinearFunction` constraints).
"""
mutable struct NonlinearBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge
    constraints::Vector
    func::MOI.VectorOfVariables
    set::ComplementsWithSetType{S}
    reformulation::AbstractComplementarityRelaxation
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{NonlinearBridge{T,S}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
    # Only `MathOptComplements.Optimizer` supports setting custom bridge arguments
    # so the default will be used when the bridge is added for instance to
    # `MOI.Bridges.LazyBridgeOptimizer`
    reformulation::AbstractComplementarityRelaxation = ScholtesRelaxation(
        zero(T),
    ),
) where {T,S}
    # Delay reformulation until `final_touch` so that per-constraint
    # `ComplementarityReformulation` attributes can override it first.
    return NonlinearBridge{T,S}([], func, set, reformulation)
end

function MOI.supports_constraint(
    ::Type{<:NonlinearBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:NonlinearBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    return NonlinearBridge{T,S}
end

function MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:NonlinearBridge},
)
    return true
end

function MOI.set(
    model::MOI.ModelLike,
    ::ComplementarityReformulation,
    bridge::NonlinearBridge,
    value::AbstractComplementarityRelaxation,
)
    bridge.reformulation = value
    # Delete old constraints so that final_touch re-reformulates
    for ci in bridge.constraints
        MOI.delete(model, ci)
    end
    empty!(bridge.constraints)
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:NonlinearBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{<:NonlinearBridge})
    return Tuple{Type,Type}[]
end

function MOI.get(
    bridge::NonlinearBridge,
    ::MOI.NumberOfConstraints{F,S},
)::Int64 where {F,S}
    return count(ci -> ci isa MOI.ConstraintIndex{F,S}, bridge.constraints)
end

function MOI.get(
    bridge::NonlinearBridge,
    ::MOI.ListOfConstraintIndices{F,S},
) where {F,S}
    return MOI.ConstraintIndex{F,S}[
        ci for ci in bridge.constraints if ci isa MOI.ConstraintIndex{F,S}
    ]
end

function MOI.delete(model::MOI.ModelLike, bridge::NonlinearBridge)
    for ci in bridge.constraints
        MOI.delete(model, ci)
    end
    return
end

MOI.Bridges.needs_final_touch(::NonlinearBridge) = true

function MOI.Bridges.final_touch(bridge::NonlinearBridge, model::MOI.ModelLike)
    if !isempty(bridge.constraints)
        return
    end
    append!(
        bridge.constraints,
        _reformulate!(model, bridge.reformulation, bridge.func, bridge.set),
    )
    return
end

# Bound helpers: extract the relevant bounds based on the set type S.
function _complementarity_bounds(
    ::Type{MOI.Nonnegatives},
    model,
    ::Type{T},
    x,
) where {T}
    return (zero(T), typemax(T))
end

function _complementarity_bounds(
    ::Type{MOI.Nonpositives},
    model,
    ::Type{T},
    x,
) where {T}
    return (typemin(T), zero(T))
end

function _complementarity_bounds(
    ::Type{MOI.Zeros},
    model,
    ::Type{T},
    x,
) where {T}
    return (zero(T), zero(T))
end

function _complementarity_bounds(
    ::Type{MOI.GreaterThan{T}},
    model,
    ::Type{T},
    x,
) where {T}
    return (MOI.Utilities.get_bounds(model, T, x)[1], typemax(T))
end

function _complementarity_bounds(
    ::Type{MOI.LessThan{T}},
    model,
    ::Type{T},
    x,
) where {T}
    return (typemin(T), MOI.Utilities.get_bounds(model, T, x)[2])
end

function _complementarity_bounds(
    ::Type{MOI.Interval{T}},
    model,
    ::Type{T},
    x,
) where {T}
    return MOI.Utilities.get_bounds(model, T, x)
end

function _reformulate!(
    model::MOI.ModelLike,
    relaxation::AbstractComplementarityRelaxation,
    fun::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
) where {S}
    n_comp = div(set.dimension, 2)
    ret = MOI.ConstraintIndex[]
    for cc in 1:n_comp
        F, x = fun.variables[cc], fun.variables[cc+n_comp]
        l, u = _complementarity_bounds(S, model, Float64, x)
        _reformulate!(ret, model, relaxation, F, x, l, u)
    end
    return ret
end

"""
    ScholtesRelaxation{T<:Real}(tau::T) <: AbstractComplementarityRelaxation

Implement the Scholtes relaxation.

For ``\\tau \\ge 0``, the complementarity constraint
```math
f(x) \\perp l \\le x \\le u
```
is reformulated as
```math
\\begin{aligned}
f(x) \\cdot (x - l) \\le \\tau       \\\\
f(x) \\cdot (x - u) \\le \\tau
\\end{aligned}
```
If ``\\tau`` is equal to 0, then the relaxation is exact.
"""
struct ScholtesRelaxation{T<:Real} <: AbstractComplementarityRelaxation
    tau::T
end

function _reformulate!(
    ret::Vector{MOI.ConstraintIndex},
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    F::MOI.VariableIndex,
    x::MOI.VariableIndex,
    l::T,
    u::T,
) where {T}
    if !isinf(l)
        f = F * (x - l)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(relaxation.tau)))
    end
    if !isinf(u)
        f = F * (x - u)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(relaxation.tau)))
    end
    return
end

"""
    FischerBurmeisterRelaxation{T}(epsilon::T) <:
    AbstractComplementarityRelaxation

Implement the Fischer-Burmeister relaxation for complementarity constraints.

For ``\\epsilon \\ge 0``, the complementarity constraint
```math
f(x) \\perp l \\le x \\le u
```
is reformulated as
```math
\\begin{aligned}
\\phi(u(x), x - l) \\le 0       \\\\
\\phi(-f(x), u - x) \\le 0
\\end{aligned}
```
where
```math
\\phi(a, b) = a + b - \\sqrt{a^2 + b^2 + \\epsilon}
```
"""
struct FischerBurmeisterRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _phi(relaxation::FischerBurmeisterRelaxation, a, b)
    eps = relaxation.epsilon
    return MOI.ScalarNonlinearFunction(
        :-,
        Any[
            1.0*a+1.0*b,
            MOI.ScalarNonlinearFunction(:sqrt, Any[(1.0*a)^2+(1.0*b)^2+eps]),
        ],
    )
end

function _reformulate!(
    ret::Vector{MOI.ConstraintIndex},
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    F::MOI.VariableIndex,
    x::MOI.VariableIndex,
    l::T,
    u::T,
) where {T}
    if !isinf(l)
        f = _phi(relaxation, F, x - l)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(0.0)))
    end
    if !isinf(u)
        f = _phi(relaxation, -1.0 * F, u - x)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(0.0)))
    end
    return
end

"""
    LiuFukushimaRelaxation <: AbstractComplementarityRelaxation

Implement the Liu-Fukushima relaxation for complementarity constraints.

For ``\\epsilon \\ge 0``, the complementarity constraint:
```math
f(x) \\perp l \\le x \\le u
```
is reformulated, if ``l > -\\infty`` as:
```math
\\begin{aligned}
f(x) \\cdot (x - l) \\le \\epsilon^2                                \\\\
(f(x) + \\epsilon) \\cdot (x - l + \\epsilon) \\ge \\epsilon^2
```
or, if ``u < \\infty``:
```math
-f(x) \\cdot (u - x) \\le \\epsilon^2                               \\\\
(-f(x) + \\epsilon) \\cdot (u - x + \\epsilon) \\ge \\epsilon^2
\\end{aligned}
```

This reformulation does not support the case where
``-\\infty < l \\le u < \\infty``.
"""
struct LiuFukushimaRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _reformulate!(
    ret::Vector{MOI.ConstraintIndex},
    model::MOI.ModelLike,
    relaxation::LiuFukushimaRelaxation,
    F::MOI.VariableIndex,
    x::MOI.VariableIndex,
    l::T,
    u::T,
) where {T}
    @assert isinf(l) || isinf(u)
    ε = relaxation.epsilon
    if !isinf(l)
        push!(ret, MOI.add_constraint(model, F * (x - l), MOI.LessThan(ε^2)))
        f = (F + ε) * (x - l + ε)
        push!(ret, MOI.add_constraint(model, f, MOI.GreaterThan(ε^2)))
    elseif !isinf(u)
        push!(ret, MOI.add_constraint(model, F * (x - u), MOI.LessThan(ε^2)))
        f = (F - ε) * (x - u - ε)
        push!(ret, MOI.add_constraint(model, f, MOI.GreaterThan(ε^2)))
    end
    return
end

"""
    KanzowSchwarzRelaxation <: AbstractComplementarityRelaxation

Implement the Kanzow-Schwarz relaxation for complementarity constraints.

For ``\\epsilon \\ge 0``, the complementarity constraint
```math
f(x) \\perp l \\le x \\le u
```
is reformulated as
```math
\\begin{aligned}
\\phi(f(x), x - l) \\le 0       \\\\
\\phi(-f(x), u - x) \\le 0
\\end{aligned}
```
where
```
ϕ(a, b) = (a - ε) * (b - ε)           if a + b > 2ε
          -((a - ε)^2 + (b - ε)^2)/2  otherwise
```
"""
struct KanzowSchwarzRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _phi(relaxation::KanzowSchwarzRelaxation, a, b)
    eps = relaxation.epsilon
    return MOI.ScalarNonlinearFunction(
        :ifelse,
        Any[
            MOI.ScalarNonlinearFunction(:>, Any[1.0*a+1.0*b, 2*eps]),
            (a-eps)*(b-eps),
            -0.5*((a-eps)^2+(b-eps)^2),
        ],
    )
end

function _reformulate!(
    ret::Vector{MOI.ConstraintIndex},
    model::MOI.ModelLike,
    relaxation::KanzowSchwarzRelaxation,
    F::MOI.VariableIndex,
    x::MOI.VariableIndex,
    l::T,
    u::T,
) where {T}
    if !isinf(l)
        f = _phi(relaxation, F, x - l)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(0.0)))
    end
    if !isinf(u)
        f = _phi(relaxation, -1.0 * F, u - x)
        push!(ret, MOI.add_constraint(model, f, MOI.LessThan(0.0)))
    end
    return
end
