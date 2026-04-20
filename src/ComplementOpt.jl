module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("utils.jl")
include("attributes.jl")
include("sets.jl")
include("vertical.jl")
include("specify_set_type_bridge.jl")
include("complements_vectorize_bridge.jl")
include("nonlinear.jl")
include("sos1.jl")
include("MOI_wrapper.jl")

end # module ComplementOpt
