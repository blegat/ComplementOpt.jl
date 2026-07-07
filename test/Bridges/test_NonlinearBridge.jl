# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestNonlinearBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestNonlinearBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_Nonnegatives_lower_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        """
        variables: x, y
        y >= 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """,
        """
        variables: x, y
        y >= 0.0
        1.0 * x * y <= 0.0
        """;
        cannot_unbridge = true,
    )
    return
end

function test_Nonpositives_upper_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        """
        variables: x, y
        y <= 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Nonpositives}(2)
        """,
        """
        variables: x, y
        y <= 0.0
        1.0 * x * y <= 0.0
        """;
        cannot_unbridge = true,
    )
    return
end

function test_GreaterThan_with_nonzero_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        """
        variables: x, y
        y >= 3.0
        [x, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        """,
        """
        variables: x, y
        y >= 3.0
        1.0 * x * y + -3.0 * x <= 0.0
        """;
        cannot_unbridge = true,
    )
    return
end

function test_Nonnegatives_with_unbounded_x1_lower_bound()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    MOI.Bridges.final_touch(model)
    target = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    x1t = MOI.add_variable(target)
    x2t, _ = MOI.add_constrained_variable(target, MOI.GreaterThan(0.0))
    MOI.add_constraint(target, x1t * (x2t - 0.0), MOI.LessThan(0.0))
    MOI.Bridges._test_structural_identical(target, inner)
    return
end

function test_Nonpositives_with_unbounded_x1_upper_bound()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
    )
    MOI.Bridges.final_touch(model)
    target = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    x1t = MOI.add_variable(target)
    x2t, _ = MOI.add_constrained_variable(target, MOI.LessThan(0.0))
    MOI.add_constraint(target, x1t * (x2t - 0.0), MOI.LessThan(0.0))
    MOI.Bridges._test_structural_identical(target, inner)
    return
end

function test_LessThan_with_nonzero_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        """
        variables: x, y
        y <= -2.0
        [x, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """,
        """
        variables: x, y
        y <= -2.0
        1.0 * x * y + 2.0 * x <= 0.0
        """;
        cannot_unbridge = true,
    )
    return
end

function test_FB_Nonnegatives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, MathOptComplements.ComplementarityReformulation(), ci, relax)
    MOI.Bridges.final_touch(model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VariableIndex,MOI.GreaterThan{Float64}}(),
    ) == 1
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.LessThan{Float64},
        }(),
    ) == 1
    return
end

function test_FB_Nonpositives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
    )
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, MathOptComplements.ComplementarityReformulation(), ci, relax)
    MOI.Bridges.final_touch(model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VariableIndex,MOI.LessThan{Float64}}(),
    ) == 1
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.LessThan{Float64},
        }(),
    ) == 1
    return
end

function test_KS_Nonnegatives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    relax = MathOptComplements.KanzowSchwarzRelaxation(1e-8)
    MOI.set(model, MathOptComplements.ComplementarityReformulation(), ci, relax)
    MOI.Bridges.final_touch(model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VariableIndex,MOI.GreaterThan{Float64}}(),
    ) == 1
    return
end

function test_KS_Nonpositives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
    )
    relax = MathOptComplements.KanzowSchwarzRelaxation(1e-8)
    MOI.set(model, MathOptComplements.ComplementarityReformulation(), ci, relax)
    MOI.Bridges.final_touch(model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VariableIndex,MOI.LessThan{Float64}}(),
    ) == 1
    return
end

# `MOI.Bridges.runtests` cannot catch an under-specified
# `added_constraint_types`: its `_general_bridge_tests` only iterates over the
# *declared* types to check consistency, so it never verifies that every
# constraint actually added by the bridge is declared. We check that here
# explicitly for each relaxation, since the produced constraint types depend on
# the relaxation (quadratic for Scholtes/Liu-Fukushima, nonlinear for
# Fischer-Burmeister/Kanzow-Schwarz).
function _test_added_constraint_types(relax)
    T = Float64
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    B = MathOptComplements.Bridges.NonlinearBridge{T}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(zero(T)))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    if relax !== nothing
        MOI.set(
            model,
            MathOptComplements.ComplementarityReformulation(),
            ci,
            relax,
        )
    end
    MOI.Bridges.final_touch(model)
    CB = MathOptComplements.Bridges.NonlinearBridge{T,MOI.Nonnegatives}
    # The bridge adds no variable, so `added_constrained_variable_types` is
    # empty.
    @test isempty(MOI.Bridges.added_constrained_variable_types(CB))
    declared = MOI.Bridges.added_constraint_types(CB)
    for (F, S) in MOI.get(inner, MOI.ListOfConstraintTypesPresent())
        # The bound on `x2` is added by the test, not by the bridge.
        if F == MOI.VariableIndex
            continue
        end
        @test (F, S) in declared
    end
    return
end

function test_added_constraint_types()
    _test_added_constraint_types(nothing)  # ScholtesRelaxation (the default)
    _test_added_constraint_types(
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    _test_added_constraint_types(
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
    )
    _test_added_constraint_types(
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    )
    return
end

function test_Zeros_equality()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        """
        variables: x, y
        [x, y] in $M.ComplementsWithSetType{MOI.Zeros}(2)
        """,
        """
        variables: x, y
        1.0 * x * y <= 0.0
        1.0 * x * y <= 0.0
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestNonlinearBridge

TestNonlinearBridge.runtests()
