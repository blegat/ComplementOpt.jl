module MathOptComplementsJuMPExt

import MathOptComplements
import JuMP
import MathOptInterface as MOI

function MathOptComplements.Bridges.add_all_bridges(
    model::JuMP.GenericModel{T},
    ::Type{U} = T,
) where {T,U}
    for BT in (
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        MathOptComplements.Bridges.ComplementsVectorizeBridge,
        MathOptComplements.Bridges.SplitIntervalBridge,
        MathOptComplements.Bridges.FlipSignBridge,
        MathOptComplements.Bridges.ToSOS1Bridge,
        MathOptComplements.Bridges.VerticalBridge,
        MathOptComplements.Bridges.NonlinearBridge,
    )
        JuMP.add_bridge(model, BT; coefficient_type = U)
    end
    return
end

end # module MathOptComplementsJuMPExt
