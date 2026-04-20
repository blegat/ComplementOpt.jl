"""
    SpecifySetTypeBridge{T} <: MOI.Bridges.Constraint.AbstractBridge

Bridge that converts a `VectorOfVariables`-in-`Complements` constraint into
per-pair `VectorOfVariables`-in-`ComplementsWithSetType{S}` constraints, where
`S` is determined by reading each slack variable's bounds.

No slack variables are created — this bridge only classifies each pair and adds
the appropriate bound on the activity variable `x₁`.

"""
mutable struct SpecifySetTypeBridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    constraints::Vector{MOI.ConstraintIndex}
    bounds::Vector{MOI.ConstraintIndex}
    func::MOI.VectorOfVariables
    set::MOI.Complements
    reformulation::Union{Nothing,AbstractComplementarityRelaxation}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SpecifySetTypeBridge{T}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::MOI.Complements,
) where {T}
    return SpecifySetTypeBridge{T}(MOI.ConstraintIndex[], MOI.ConstraintIndex[], func, set, nothing)
end

function MOI.supports_constraint(
    ::Type{<:SpecifySetTypeBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Complements},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{SpecifySetTypeBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Complements},
) where {T}
    return SpecifySetTypeBridge{T}
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

function MOI.Bridges.final_touch(
    bridge::SpecifySetTypeBridge{T},
    model::MOI.ModelLike,
) where {T}
    if !isempty(bridge.constraints)
        return
    end
    n_comp = div(bridge.set.dimension, 2)
    for cc in 1:n_comp
        x1 = bridge.func.variables[cc]
        x2 = bridge.func.variables[cc+n_comp]
        ci = _specify_set_type_pair!(model, T, x1, x2, bridge.bounds)
        push!(bridge.constraints, ci)
    end
    if bridge.reformulation !== nothing
        for ci in bridge.constraints
            MOI.set(model, ComplementarityReformulation(), ci, bridge.reformulation)
        end
    end
    return
end

function _specify_set_type_pair!(model, ::Type{T}, x1, x2, bounds) where {T}
    lb2, ub2 = MOIU.get_bounds(model, T, x2)
    if !isinf(lb2) && isinf(ub2)
        return _specify_lower_bound!(model, T, x1, x2, lb2, bounds)
    elseif isinf(lb2) && !isinf(ub2)
        return _specify_upper_bound!(model, T, x1, x2, ub2, bounds)
    elseif isfinite(lb2) && isfinite(ub2)
        return _specify_range!(model, T, x1, x2, lb2, ub2)
    else
        # Both infinite: x1 must be zero
        push!(bounds, MOI.add_constraint(model, one(T) * x1, MOI.EqualTo(zero(T))))
        return MOI.add_constraint(
            model,
            MOI.VectorOfVariables([x1, x2]),
            ComplementsWithSetType{MOI.Zeros}(2),
        )
    end
end

function _specify_lower_bound!(model, ::Type{T}, x1, x2, lb2, bounds) where {T}
    lb1, _ = MOIU.get_bounds(model, T, x1)
    if isinf(lb1)
        push!(bounds, MOI.add_constraint(model, x1, MOI.GreaterThan(zero(T))))
    end
    S = iszero(lb2) ? MOI.Nonnegatives : MOI.GreaterThan{T}
    return MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        ComplementsWithSetType{S}(2),
    )
end

function _specify_upper_bound!(model, ::Type{T}, x1, x2, ub2, bounds) where {T}
    _, ub1 = MOIU.get_bounds(model, T, x1)
    if isinf(ub1)
        push!(bounds, MOI.add_constraint(model, x1, MOI.LessThan(zero(T))))
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

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:SpecifySetTypeBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{SpecifySetTypeBridge{T}},
) where {T}
    return Tuple{Type,Type}[
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.Nonnegatives}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.Nonpositives}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.Zeros}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.GreaterThan{T}}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.LessThan{T}}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.Interval{T}}),
        (MOI.VariableIndex, MOI.GreaterThan{T}),
        (MOI.VariableIndex, MOI.LessThan{T}),
        (MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}),
    ]
end

function MOI.get(::SpecifySetTypeBridge, ::MOI.NumberOfVariables)::Int64
    return 0
end

function MOI.get(::SpecifySetTypeBridge, ::MOI.ListOfVariableIndices)
    return MOI.VariableIndex[]
end

function MOI.get(
    bridge::SpecifySetTypeBridge,
    ::MOI.NumberOfConstraints{F,S},
)::Int64 where {F,S}
    all_cis = [bridge.constraints; bridge.bounds]
    return count(ci -> ci isa MOI.ConstraintIndex{F,S}, all_cis)
end

function MOI.get(
    bridge::SpecifySetTypeBridge,
    ::MOI.ListOfConstraintIndices{F,S},
) where {F,S}
    all_cis = [bridge.constraints; bridge.bounds]
    return MOI.ConstraintIndex{F,S}[
        ci for ci in all_cis if ci isa MOI.ConstraintIndex{F,S}
    ]
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::SpecifySetTypeBridge,
)
    return bridge.func
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::SpecifySetTypeBridge,
)
    return bridge.set
end

function MOI.delete(model::MOI.ModelLike, bridge::SpecifySetTypeBridge)
    for ci in bridge.constraints
        MOI.delete(model, ci)
    end
    for ci in bridge.bounds
        MOI.delete(model, ci)
    end
    return
end
