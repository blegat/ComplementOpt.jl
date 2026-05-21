# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestVerticalBridge

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

function test_Unbounded_RHS_with_constant_in_LHS_Complements()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.VerticalBridge,
        model -> begin
            x = MOI.add_variable(model)
            z = MOI.add_variable(model)
            y = MOI.add_variable(model)
            w, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            f = MOI.VectorAffineFunction{Float64}(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, z)),
                    MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(1.0, y)),
                    MOI.VectorAffineTerm(4, MOI.ScalarAffineTerm(1.0, w)),
                ],
                [3.0, 0.0, 0.0, 0.0],
            )
            MOI.add_constraint(model, f, MOI.Complements(4))
        end,
        model -> begin
            x = MOI.add_variable(model)
            z = MOI.add_variable(model)
            y = MOI.add_variable(model)
            w, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            # 2x + 3 = 0  →  2x in EqualTo(-3.0)  (y is unbounded)
            MOI.add_constraint(model, 2.0 * x, MOI.EqualTo(-3.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([z, w]),
                MOI.Complements(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_Expression_LHS_with_constant_ComplementsWithSetType()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.VerticalBridge,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            f = MOI.VectorAffineFunction{Float64}(
                [
                    MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x)),
                    MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, y)),
                ],
                [3.0, 0.0],
            )
            MOI.add_constraint(
                model,
                f,
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end,
        model -> begin
            x = MOI.add_variable(model)
            y, _ =
                MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            x1 = MOI.add_variable(model)
            # 2x + 3 = x1  →  2x - x1 in EqualTo(-3.0)
            MOI.add_constraint(model, 2.0 * x - 1.0 * x1, MOI.EqualTo(-3.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x1, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # TestVerticalBridge

TestVerticalBridge.runtests()
