module Bridges

import MathOptInterface as MOI
const MOIU = MOI.Utilities
using ..ComplementOpt:
    ComplementsWithSetType,
    AbstractComplementarityRelaxation,
    ComplementarityReformulation,
    _remove_bounds!

include("vertical.jl")
include("specify_set_type_bridge.jl")
include("complements_vectorize_bridge.jl")
include("split_interval_bridge.jl")
include("flip_sign_bridge.jl")
include("nonlinear.jl")
include("to_sos1_bridge.jl")
include("sos1.jl")

"""
    add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64)

Add all `ComplementOpt` bridges to `model`. The model is typically a
[`MOI.Bridges.LazyBridgeOptimizer`](@ref) so that the bridge graph is
extended with the bridges needed to reformulate
[`ComplementOpt.ComplementsWithSetType`](@ref) and [`MOI.Complements`](@ref)
constraints.

When used with a `LazyBridgeOptimizer`, the [`NonlinearBridge`](@ref) uses
the default [`ScholtesRelaxation`](@ref) because the
[`ComplementOpt.DefaultComplementarityReformulation`](@ref) optimizer
attribute is only supported by [`ComplementOpt.Optimizer`](@ref).
"""
function add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64) where {T}
    MOI.Bridges.add_bridge(model, SpecifySetTypeBridge{T})
    MOI.Bridges.add_bridge(model, ComplementsVectorizeBridge{T})
    MOI.Bridges.add_bridge(model, SplitIntervalBridge{T})
    MOI.Bridges.add_bridge(model, FlipSignBridge{T})
    MOI.Bridges.add_bridge(model, ToSOS1Bridge{T})
    MOI.Bridges.add_bridge(model, VerticalBridge{T})
    MOI.Bridges.add_bridge(model, NonlinearBridge{T})
    return
end

end # module Bridges
