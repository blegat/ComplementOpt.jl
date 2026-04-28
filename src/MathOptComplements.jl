module MathOptComplements

import MathOptInterface as MOI
const MOIU = MOI.Utilities

include("utils.jl")
include("attributes.jl")
include("sets.jl")
include("Bridges/Bridges.jl")

using .Bridges:
    ScholtesRelaxation,
    FischerBurmeisterRelaxation,
    LiuFukushimaRelaxation,
    KanzowSchwarzRelaxation

include("MOI_wrapper.jl")

end # module MathOptComplements
