# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestComplementsVectorizeBridge

using Test
using JuMP
using MathOptComplements

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_GreaterThan_Nonnegatives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
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
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
            f = MOI.VectorAffineFunction{Float64}(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, y)),
                ],
                [0.0, -3.0],
            )
            MOI.add_constraint(
                model,
                f,
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_LessThan_Nonpositives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
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
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
            f = MOI.VectorAffineFunction{Float64}(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, y)),
                ],
                [0.0, 2.0],
            )
            MOI.add_constraint(
                model,
                f,
                MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_EqualTo_Zeros()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.EqualTo(2.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.EqualTo{Float64}}(
                    2,
                ),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.EqualTo(2.0))
            f = MOI.Utilities.operate(vcat, Float64, 1.0 * x, 1.0 * y - 2.0)
            MOI.add_constraint(
                model,
                f,
                MathOptComplements.ComplementsWithSetType{MOI.Zeros}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # TestComplementsVectorizeBridge

TestComplementsVectorizeBridge.runtests()
