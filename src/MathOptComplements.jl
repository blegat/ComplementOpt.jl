# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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
