
function _remove_bounds!(model::MOI.ModelLike, x::MOI.VariableIndex)
    for cidx in [
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}}(x.value),
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x.value),
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x.value),
    ]
        if MOI.is_valid(model, cidx)
            MOI.delete(model, cidx)
        end
    end
    return
end
