struct Reformulation <: MOI.AbstractOptimizerAttribute end

struct Optimizer{O<:MOI.ModelLike} <: MOI.AbstractOptimizer
    optimizer::O
    reformulation
    function Optimizer(optimizer::MOI.ModelLike)
        return new{typeof(optimizer)}(
            optimizer,
            NonlinearReformulation(),
        )
    end
end

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.optimizer)
MOI.empty!(model::Optimizer) = MOI.empty!(model.optimizer)

function MOI.supports(
    model::Optimizer,
    attr::Union{
        MOI.AbstractModelAttribute,
        MOI.AbstractOptimizerAttribute,
    },
)
    return MOI.supports(model.optimizer, attr)
end

function MOI.get(
    model::Optimizer,
    attr::Union{
        MOI.AbstractModelAttribute,
        MOI.AbstractOptimizerAttribute,
    },
)
    return MOI.get(model.optimizer, attr)
end

function MOI.set(
    model::Optimizer,
    attr::Union{
        MOI.AbstractModelAttribute,
        MOI.AbstractOptimizerAttribute,
    },
    value,
)
    return MOI.set(model.optimizer, attr, value)
end

function MOI.supports_constraint(
    model::Optimizer,
    ::Type{F},
    ::Type{S},
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return MOI.supports_constraint(model.optimizer, F, S)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    tmp = MOI.Utilities.UniversalFallback(
        MOI.Utilities.Model{Float64}()
    )
    tmp_index_map = MOI.copy_to(tmp, src)
    reformulate_to_vertical!(tmp)
    reformulate_as_nonlinear_program!(tmp)
    # TODO combine with `tmp_index_map`
    return MOI.copy_to(dest.optimizer, tmp)
end
