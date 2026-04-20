module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("utils.jl")
include("attributes.jl")
include("sets.jl")
include("Bridges/Bridges.jl")

# Re-export bridge types and relaxation types at the ComplementOpt level
using .Bridges:
    VerticalBridge,
    SpecifySetTypeBridge,
    ComplementsVectorizeBridge,
    NonlinearBridge,
    ScholtesRelaxation,
    FischerBurmeisterRelaxation,
    LiuFukushimaRelaxation,
    KanzowSchwarzRelaxation,
    SOS1Relaxation

include("MOI_wrapper.jl")

end # module ComplementOpt
