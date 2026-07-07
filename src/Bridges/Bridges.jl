# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Bridges

import MathOptInterface as MOI

using ..MathOptComplements:
    ComplementsWithSetType,
    AbstractComplementarityRelaxation,
    ComplementarityReformulation

include("VerticalBridge.jl")
include("SpecifySetTypeBridge.jl")
include("ComplementsVectorizeBridge.jl")
include("SplitIntervalBridge.jl")
include("FlipSignBridge.jl")
include("NonlinearBridge.jl")
include("ToSOS1Bridge.jl")
include("ToZerosBridge.jl")

const _ALL_BRIDGE_TYPES = Any[
    SpecifySetTypeBridge,
    ComplementsVectorizeBridge,
    SplitIntervalBridge,
    FlipSignBridge,
    ToSOS1Bridge,
    ToZerosBridge,
    VerticalBridge,
    NonlinearBridge,
]

"""
    add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64)

Add all `MathOptComplements` bridges to `model`. The model is typically a
[`MOI.Bridges.LazyBridgeOptimizer`](@ref) so that the bridge graph is
extended with the bridges needed to reformulate
[`MathOptComplements.ComplementsWithSetType`](@ref) and [`MOI.Complements`](@ref)
constraints.

When used with a `LazyBridgeOptimizer`, the [`NonlinearBridge`](@ref) uses
the default [`ScholtesRelaxation`](@ref) because the
[`MathOptComplements.DefaultComplementarityReformulation`](@ref) optimizer
attribute is only supported by [`MathOptComplements.Optimizer`](@ref).
"""
function add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64) where {T}
    for bridge_type in _ALL_BRIDGE_TYPES
        MOI.Bridges.add_bridge(model, bridge_type{T})
    end
    return
end

end # module Bridges
