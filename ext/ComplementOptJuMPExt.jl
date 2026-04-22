module ComplementOptJuMPExt

import ComplementOpt
import JuMP
import MathOptInterface as MOI

function ComplementOpt.add_all_bridges(
    model::JuMP.GenericModel{T},
    ::Type{U} = T,
) where {T,U}
    JuMP.add_bridge(model, ComplementOpt.Bridges.SpecifySetTypeBridge; coefficient_type = U)
    JuMP.add_bridge(
        model,
        ComplementOpt.Bridges.ComplementsVectorizeBridge;
        coefficient_type = U,
    )
    JuMP.add_bridge(model, ComplementOpt.Bridges.SplitIntervalBridge; coefficient_type = U)
    JuMP.add_bridge(model, ComplementOpt.Bridges.FlipSignBridge; coefficient_type = U)
    JuMP.add_bridge(model, ComplementOpt.Bridges.ToSOS1Bridge; coefficient_type = U)
    JuMP.add_bridge(
        model,
        ComplementOpt.Bridges.VerticalBridge{MOI.Complements};
        coefficient_type = U,
    )
    JuMP.add_bridge(
        model,
        ComplementOpt.Bridges.VerticalBridge{ComplementOpt.ComplementsWithSetType};
        coefficient_type = U,
    )
    JuMP.add_bridge(model, ComplementOpt.Bridges.NonlinearBridge; coefficient_type = U)
    return
end

end # module ComplementOptJuMPExt
