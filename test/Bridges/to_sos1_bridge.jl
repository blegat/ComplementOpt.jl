# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Test
using JuMP
using MathOptComplements

@testset "ToSOS1Bridge" begin
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ToSOS1Bridge,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
            )
        end,
        model -> begin
            x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constraint(model, MOI.VectorOfVariables([x, y]), MOI.SOS1([1.0, 2.0]))
        end;
        cannot_unbridge = true,
    )
end
