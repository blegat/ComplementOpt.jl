# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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
