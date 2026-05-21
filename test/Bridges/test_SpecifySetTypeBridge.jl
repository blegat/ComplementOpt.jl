# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestSpecifySetTypeBridge

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

function test_lower_bound_Nonnegatives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_lower_bound_GreaterThan()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{
                    MOI.GreaterThan{Float64},
                }(
                    2,
                ),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_upper_bound_Nonpositives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_upper_bound_LessThan()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(
                    2,
                ),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_range_Interval()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(
                    2,
                ),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_Free_variable_Zeros()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        model -> begin
            x = MOI.add_variable(model)
            y = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.Complements(2),
            )
        end,
        model -> begin
            x = MOI.add_variable(model)
            y = MOI.add_variable(model)
            # x1 must be zero
            MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Zeros}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # TestSpecifySetTypeBridge

TestSpecifySetTypeBridge.runtests()
