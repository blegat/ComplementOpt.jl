# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestNonlinearBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

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
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x2, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, x2]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x2, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            # x1 * (x2 - 0) <= 0
            MOI.add_constraint(model, x1 * (x2 - 0.0), MOI.LessThan(0.0))
        end;
        cannot_unbridge = true,
    )
    return
end

function test_Nonpositives_upper_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, x2]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
            )
        end,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            # x1 * (x2 - 0) <= 0
            MOI.add_constraint(model, x1 * (x2 - 0.0), MOI.LessThan(0.0))
        end;
        cannot_unbridge = true,
    )
    return
end

function test_GreaterThan_with_nonzero_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x2, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, x2]),
                MathOptComplements.ComplementsWithSetType{
                    MOI.GreaterThan{Float64},
                }(
                    2,
                ),
            )
        end,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x2, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
            # x1 * (x2 - 3) <= 0
            MOI.add_constraint(model, x1 * (x2 - 3.0), MOI.LessThan(0.0))
        end;
        cannot_unbridge = true,
    )
    return
end

function test_Nonnegatives_with_unbounded_x1_lower_bound()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
    )
    MOI.Bridges.final_touch(model)
    target = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    x1t, _ = MOI.add_constrained_variable(target, MOI.GreaterThan(0.0))
    x2t, _ = MOI.add_constrained_variable(target, MOI.GreaterThan(0.0))
    MOI.add_constraint(target, x1t * (x2t - 0.0), MOI.LessThan(0.0))
    MOI.Bridges._test_structural_identical(target, inner)
    return
end

function test_Nonpositives_with_unbounded_x1_upper_bound()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
    x1 = MOI.add_variable(model)
    x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x1, x2]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
    )
    MOI.Bridges.final_touch(model)
    target = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    x1t, _ = MOI.add_constrained_variable(target, MOI.LessThan(0.0))
    x2t, _ = MOI.add_constrained_variable(target, MOI.LessThan(0.0))
    MOI.add_constraint(target, x1t * (x2t - 0.0), MOI.LessThan(0.0))
    MOI.Bridges._test_structural_identical(target, inner)
    return
end

function test_LessThan_with_nonzero_bound()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, x2]),
                MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(
                    2,
                ),
            )
        end,
        model -> begin
            x1, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
            # x1 * (x2 - (-2)) <= 0
            MOI.add_constraint(model, x1 * (x2 + 2.0), MOI.LessThan(0.0))
        end;
        cannot_unbridge = true,
    )
    return
end

function test_FB_Nonnegatives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
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
    ) == 2
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
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
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
    ) == 2
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.GreaterThan{Float64},
        }(),
    ) == 1
    return
end

function test_KS_Nonnegatives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
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
    ) == 2
    return
end

function test_KS_Nonpositives_with_unbounded_x1()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.NonlinearBridge{Float64},
    }(
        inner,
    )
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
    ) == 2
    return
end

function test_Zeros_equality()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.NonlinearBridge,
        model -> begin
            x1 = MOI.add_variable(model)
            x2 = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, x2]),
                MathOptComplements.ComplementsWithSetType{MOI.Zeros}(2),
            )
        end,
        model -> begin
            x1 = MOI.add_variable(model)
            x2 = MOI.add_variable(model)
            # Range relaxation with lb=0, ub=0: two constraints
            MOI.add_constraint(model, x1 * (x2 - 0.0), MOI.LessThan(0.0))
            MOI.add_constraint(model, x1 * (x2 - 0.0), MOI.LessThan(0.0))
        end;
        cannot_unbridge = true,
    )
    return
end

end  # TestNonlinearBridge

TestNonlinearBridge.runtests()
