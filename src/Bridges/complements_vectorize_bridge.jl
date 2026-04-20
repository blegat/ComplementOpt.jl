"""
    ComplementsVectorizeBridge{T,S,SV} <: MOI.Bridges.Constraint.AbstractBridge

Bridge that shifts the slack variable in a complementarity pair to remove
the bound constant, converting scalar-set-typed complements to vector-set-typed.

For example, `VectorOfVariables([x1, x2])`-in-`ComplementsWithSetType{GreaterThan{T}}`
becomes `VectorAffineFunction([x1, x2 - lb])`-in-`ComplementsWithSetType{Nonnegatives}`.

Inspired by `MOI.Bridges.Constraint.VectorizeBridge`.

The type parameters are:
- `T`: coefficient type
- `S`: input scalar set type (e.g., `GreaterThan{T}`)
- `SV`: output vector set type (e.g., `Nonnegatives`)

"""
struct ComplementsVectorizeBridge{T,S,SV} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{
        MOI.VectorAffineFunction{T},
        ComplementsWithSetType{SV},
    }
    set_constant::T
end

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
    ::Type{ComplementsVectorizeBridge{T,S,SV}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
) where {T,S,SV}
    @assert set.dimension == 2
    x1 = func.variables[1]
    x2 = func.variables[2]
    c = _set_constant(T, model, set, x2)
    vec_f = MOI.VectorAffineFunction{T}(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(one(T), x1)),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(one(T), x2)),
        ],
        T[zero(T), -c],
    )
    ci = MOI.add_constraint(model, vec_f, ComplementsWithSetType{SV}(2))
    return ComplementsVectorizeBridge{T,S,SV}(ci, c)
end

function MOI.supports_constraint(
    ::Type{<:ComplementsVectorizeBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ComplementsWithSetType{<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ComplementsVectorizeBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}
    SV = _vector_set_type(S)
    return ComplementsVectorizeBridge{T,S,SV}
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

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:ComplementsVectorizeBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{ComplementsVectorizeBridge{T,S,SV}},
) where {T,S,SV}
    return Tuple{Type,Type}[
        (MOI.VectorAffineFunction{T}, ComplementsWithSetType{SV}),
    ]
end

function MOI.get(::ComplementsVectorizeBridge, ::MOI.NumberOfVariables)::Int64
    return 0
end

function MOI.get(::ComplementsVectorizeBridge, ::MOI.ListOfVariableIndices)
    return MOI.VariableIndex[]
end

function MOI.get(
    ::ComplementsVectorizeBridge{T,S,SV},
    ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T},ComplementsWithSetType{SV}},
)::Int64 where {T,S,SV}
    return 1
end

function MOI.get(
    bridge::ComplementsVectorizeBridge{T,S,SV},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorAffineFunction{T},
        ComplementsWithSetType{SV},
    },
) where {T,S,SV}
    return [bridge.constraint]
end

function MOI.get(
    model::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::ComplementsVectorizeBridge,
)
    f = MOI.get(model, MOI.ConstraintFunction(), bridge.constraint)
    terms = f.terms
    vars = [t.scalar_term.variable for t in terms]
    return MOI.VectorOfVariables(vars)
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    ::ComplementsVectorizeBridge{T,S},
) where {T,S}
    return ComplementsWithSetType{S}(2)
end

function MOI.delete(model::MOI.ModelLike, bridge::ComplementsVectorizeBridge)
    MOI.delete(model, bridge.constraint)
    return
end
