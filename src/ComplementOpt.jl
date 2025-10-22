module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("vertical.jl")
include("nonlinear.jl")
include("MOI_wrapper.jl")

end # module ComplementOpt
