# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Integration tests for the bridge selection performed by a
# `MOI.Bridges.LazyBridgeOptimizer` set up with `add_all_bridges` (as opposed to
# `MathOptComplements.Optimizer`, which selects the bridge explicitly). The
# `LazyBridgeOptimizer` picks a bridge based on the constraint types each bridge
# declares it adds (`added_constraint_types`) and whether they are supported by
# the solver. In particular, for a MIP solver (which supports neither quadratic
# nor nonlinear constraints), the `NonlinearBridge` must *not* be selected.

module TestLazyBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const MOC = MathOptComplements

# A MIP solver: it supports integrality and linear constraints (including
# `SOS1` only through the `SOS1ToMILP` bridge), but neither quadratic nor
# nonlinear constraints.
MOI.Utilities.@model(
    MIPModel,
    (MOI.ZeroOne, MOI.Integer),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
    (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,),
)

function MOI.supports_constraint(
    ::MIPModel,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.AbstractScalarSet},
)
    return false
end

# An NLP solver: it supports quadratic and nonlinear constraints but not `SOS1`.
MOI.Utilities.@model(
    NLPModel,
    (MOI.ZeroOne, MOI.Integer),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
    (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros),
    (),
    (),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,),
)

function MOI.supports_constraint(
    ::NLPModel,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.AbstractScalarSet},
)
    return true
end

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function _lazy(inner)
    lazy = MOI.Bridges.full_bridge_optimizer(inner, Float64)
    MOC.Bridges.add_all_bridges(lazy, Float64)
    return lazy
end

const _SET_TYPES = (
    MOI.Nonnegatives,
    MOI.Nonpositives,
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.Interval{Float64},
    MOI.Zeros,
)

function test_mip_does_not_select_nonlinear_bridge()
    lazy = _lazy(MIPModel{Float64}())
    # The `MOI.Complements` constraint is bridgeable: this requires the whole
    # SOS1/MILP chain, including `ToZerosBridge` for the `Zeros` case produced
    # by `SpecifySetTypeBridge`.
    @test MOI.supports_constraint(lazy, MOI.VectorOfVariables, MOI.Complements)
    for Sset in _SET_TYPES
        S = MOC.ComplementsWithSetType{Sset}
        @test MOI.supports_constraint(lazy, MOI.VectorOfVariables, S)
        bt = MOI.Bridges.bridge_type(lazy, MOI.VectorOfVariables, S)
        @test !(bt <: MOC.Bridges.NonlinearBridge)
    end
    # The two terminal bridges of the MILP path.
    @test MOI.Bridges.bridge_type(
        lazy,
        MOI.VectorOfVariables,
        MOC.ComplementsWithSetType{MOI.Nonnegatives},
    ) <: MOC.Bridges.ToSOS1Bridge
    @test MOI.Bridges.bridge_type(
        lazy,
        MOI.VectorOfVariables,
        MOC.ComplementsWithSetType{MOI.Zeros},
    ) <: MOC.Bridges.ToZerosBridge
    return
end

function test_nlp_selects_nonlinear_bridge()
    lazy = _lazy(NLPModel{Float64}())
    # When the solver supports nonlinear constraints, `NonlinearBridge` is
    # selected (it is cheaper than the SOS1/MILP chain since the solver does
    # not support `SOS1` natively).
    for Sset in (MOI.Nonnegatives, MOI.Interval{Float64})
        bt = MOI.Bridges.bridge_type(
            lazy,
            MOI.VectorOfVariables,
            MOC.ComplementsWithSetType{Sset},
        )
        @test bt <: MOC.Bridges.NonlinearBridge
    end
    return
end

end  # TestLazyBridge

TestLazyBridge.runtests()
