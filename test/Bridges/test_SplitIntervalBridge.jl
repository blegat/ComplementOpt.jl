# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestSplitIntervalBridge

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

function test_VectorOfVariables_input()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
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
        end,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
            xp, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            # x == xp + xn  →  x - xp - xn == 0
            MOI.add_constraint(
                model,
                1.0 * x - 1.0 * xp - 1.0 * xn,
                MOI.EqualTo(0.0),
            )
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([xp, y]),
                MathOptComplements.ComplementsWithSetType{
                    MOI.GreaterThan{Float64},
                }(
                    2,
                ),
            )
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([xn, y]),
                MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(
                    2,
                ),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_VectorAffineFunction_input()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
            f = MOI.VectorAffineFunction{Float64}(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, y)),
                ],
                [1.0, 0.0],
            )
            MOI.add_constraint(
                model,
                f,
                MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(
                    2,
                ),
            )
        end,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
            xp, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            # 2x + 1 == xp + xn  →  2x - xp - xn == -1
            MOI.add_constraint(
                model,
                2.0 * x - 1.0 * xp - 1.0 * xn,
                MOI.EqualTo(-1.0),
            )
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([xp, y]),
                MathOptComplements.ComplementsWithSetType{
                    MOI.GreaterThan{Float64},
                }(
                    2,
                ),
            )
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([xn, y]),
                MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(
                    2,
                ),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # TestSplitIntervalBridge

TestSplitIntervalBridge.runtests()
