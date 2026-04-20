"""
    ComplementsVectorizeBridge{T} <: MOI.Bridges.Constraint.AbstractBridge

Bridge that shifts the slack variable in a complementarity pair to remove
the bound constant, converting scalar-set-typed complements to vector-set-typed.

For example, `VectorOfVariables([x1, x2])`-in-`ComplementsWithSetType{GreaterThan{T}}`
becomes `VectorAffineFunction([x1, x2 - lb])`-in-`ComplementsWithSetType{Nonnegatives}`.

Inspired by `MOI.Bridges.Constraint.VectorizeBridge`.

"""
struct ComplementsVectorizeBridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex
    set_constant::T
end

const _VECTORIZE_SET_MAP = Dict(
    MOI.GreaterThan => MOI.Nonnegatives,
    MOI.LessThan => MOI.Nonpositives,
    MOI.EqualTo => MOI.Zeros,
)

function _vector_set_type(::Type{<:MOI.GreaterThan})
    return MOI.Nonnegatives
end
function _vector_set_type(::Type{<:MOI.LessThan})
    return MOI.Nonpositives
end
function _vector_set_type(::Type{<:MOI.EqualTo})
    return MOI.Zeros
end

function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.GreaterThan},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[1]
end
function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.LessThan},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[2]
end
function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.EqualTo},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[1]
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ComplementsVectorizeBridge{T}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
) where {T,S}
    @assert set.dimension == 2
    x1 = func.variables[1]
    x2 = func.variables[2]
    c = _set_constant(T, model, set, x2)
    # Build VectorAffineFunction: [1*x1 + 0, 1*x2 - c]
    vec_f = MOI.VectorAffineFunction{T}(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(one(T), x1)),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(one(T), x2)),
        ],
        T[zero(T), -c],
    )
    S_vec = _vector_set_type(S)
    ci = MOI.add_constraint(model, vec_f, ComplementsWithSetType{S_vec}(2))
    return ComplementsVectorizeBridge{T}(ci, c)
end

MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:ComplementsVectorizeBridge},
) = true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::ComplementsVectorizeBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.constraint, value)
    return
end
