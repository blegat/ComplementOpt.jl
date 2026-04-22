module ComplementOpt

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
    KanzowSchwarzRelaxation,
    SOS1Relaxation

include("MOI_wrapper.jl")

"""
    add_all_bridges(model, ::Type{T} = Float64)

Add all `ComplementOpt` bridges to `model`.

See [`ComplementOpt.Bridges.add_all_bridges`](@ref) for details.
"""
add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64) where {T} =
    Bridges.add_all_bridges(model, T)

end # module ComplementOpt
