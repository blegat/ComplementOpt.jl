module ComplementOptJuMPExt

import ComplementOpt
import JuMP
import MathOptInterface as MOI

function ComplementOpt.Bridges.add_all_bridges(
    model::JuMP.GenericModel{T},
    ::Type{U} = T,
) where {T,U}
    for BT in (
        ComplementOpt.Bridges.SpecifySetTypeBridge,
        ComplementOpt.Bridges.ComplementsVectorizeBridge,
        ComplementOpt.Bridges.SplitIntervalBridge,
        ComplementOpt.Bridges.FlipSignBridge,
        ComplementOpt.Bridges.ToSOS1Bridge,
        ComplementOpt.Bridges.VerticalBridge,
        ComplementOpt.Bridges.NonlinearBridge,
    )
        JuMP.add_bridge(model, BT; coefficient_type = U)
    end
    return
end

end # module ComplementOptJuMPExt
