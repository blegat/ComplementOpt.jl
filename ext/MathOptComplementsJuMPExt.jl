# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module MathOptComplementsJuMPExt

import MathOptComplements
import JuMP
import MathOptInterface as MOI

function MathOptComplements.Bridges.add_all_bridges(
    model::JuMP.GenericModel{T},
    ::Type{U} = T,
) where {T,U}
    for bridge_type in MathOptComplements.Bridges._ALL_BRIDGE_TYPES
        JuMP.add_bridge(model, bridge_type; coefficient_type = U)
    end
    return
end

end # module MathOptComplementsJuMPExt
