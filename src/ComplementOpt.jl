module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("utils.jl")
include("vertical.jl")
include("nonlinear.jl")
include("MOI_wrapper.jl")

end # module ComplementOpt
