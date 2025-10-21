module ComplementOpt

using JuMP
const MOIU = MOI.Utilities

include("vertical.jl")
include("nonlinear.jl")

end # module ComplementOpt
