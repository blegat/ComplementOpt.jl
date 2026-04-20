module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("utils.jl")
include("attributes.jl")
include("vertical.jl")
include("nonlinear.jl")
include("sos1.jl")
include("MOI_wrapper.jl")

end # module ComplementOpt
