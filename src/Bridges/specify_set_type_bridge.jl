"""
    SpecifySetTypeBridge <: MOI.Bridges.Constraint.AbstractBridge

Bridge that converts a `VectorOfVariables`-in-`Complements` constraint into
per-pair `VectorOfVariables`-in-`ComplementsWithSetType{S}` constraints, where
`S` is determined by reading each slack variable's bounds.

No slack variables are created — this bridge only classifies each pair and adds
the appropriate bound on the activity variable `x₁`.

"""
mutable struct SpecifySetTypeBridge <: MOI.Bridges.Constraint.AbstractBridge
    constraints::Vector{MOI.ConstraintIndex}
    func::MOI.VectorOfVariables
    set::MOI.Complements
    reformulation::Union{Nothing,AbstractComplementarityRelaxation}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SpecifySetTypeBridge},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::MOI.Complements,
)
    return SpecifySetTypeBridge(MOI.ConstraintIndex[], func, set, nothing)
end

MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:SpecifySetTypeBridge},
) = true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::SpecifySetTypeBridge,
    value::AbstractComplementarityRelaxation,
)
    bridge.reformulation = value
    for ci in bridge.constraints
        MOI.set(model, attr, ci, value)
    end
    return
end

MOI.Bridges.needs_final_touch(::SpecifySetTypeBridge) = true

function MOI.Bridges.final_touch(bridge::SpecifySetTypeBridge, model::MOI.ModelLike)
    if !isempty(bridge.constraints)
        return
    end
    n_comp = div(bridge.set.dimension, 2)
    for cc = 1:n_comp
        x1 = bridge.func.variables[cc]
        x2 = bridge.func.variables[cc+n_comp]
        ci = _specify_set_type_pair!(model, Float64, x1, x2)
        push!(bridge.constraints, ci)
    end
    if bridge.reformulation !== nothing
        for ci in bridge.constraints
            MOI.set(model, ComplementarityReformulation(), ci, bridge.reformulation)
        end
    end
    return
end

function _specify_set_type_pair!(model, ::Type{T}, x1, x2) where {T}
    lb2, ub2 = MOIU.get_bounds(model, T, x2)
    if !isinf(lb2) && isinf(ub2)
        return _specify_lower_bound!(model, T, x1, x2, lb2)
    elseif isinf(lb2) && !isinf(ub2)
        return _specify_upper_bound!(model, T, x1, x2, ub2)
    elseif isfinite(lb2) && isfinite(ub2)
        return _specify_range!(model, T, x1, x2, lb2, ub2)
    else
        # Both infinite: x1 must be zero
        MOI.add_constraint(model, 1.0 * x1, MOI.EqualTo(zero(T)))
        return MOI.add_constraint(
            model,
            MOI.VectorOfVariables([x1, x2]),
            ComplementsWithSetType{MOI.Zeros}(2),
        )
    end
end

function _specify_lower_bound!(model, ::Type{T}, x1, x2, lb2) where {T}
    lb1, _ = MOIU.get_bounds(model, T, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(zero(T)))
    end
    S = iszero(lb2) ? MOI.Nonnegatives : MOI.GreaterThan{T}
    return MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        ComplementsWithSetType{S}(2),
    )
end

function _specify_upper_bound!(model, ::Type{T}, x1, x2, ub2) where {T}
    _, ub1 = MOIU.get_bounds(model, T, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.LessThan(zero(T)))
    end
    S = iszero(ub2) ? MOI.Nonpositives : MOI.LessThan{T}
    return MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        ComplementsWithSetType{S}(2),
    )
end

function _specify_range!(model, ::Type{T}, x1, x2, lb2, ub2) where {T}
    return MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        ComplementsWithSetType{MOI.Interval{T}}(2),
    )
end
