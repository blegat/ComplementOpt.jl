module Bridges

using JuMP
const MOIU = MOI.Utilities
using ..ComplementOpt:
    ComplementsWithSetType,
    AbstractComplementarityRelaxation,
    ComplementarityReformulation,
    _remove_bounds!

include("vertical.jl")
include("specify_set_type_bridge.jl")
include("complements_vectorize_bridge.jl")
include("nonlinear.jl")
include("sos1.jl")

end # module Bridges
