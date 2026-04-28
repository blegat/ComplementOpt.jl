using Test
using JuMP
using MathOptComplements

@testset "NonlinearBridge" begin
    @testset "Nonnegatives (lower bound)" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.NonlinearBridge,
            model -> begin
                x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x1, x2]),
                    MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end,
            model -> begin
                x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                # x1 * (x2 - 0) <= 0
                MOI.add_constraint(model, x1 * (x2 - 0.0), MOI.LessThan(0.0))
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Nonpositives (upper bound)" begin
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
    end

    @testset "GreaterThan with non-zero bound" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.NonlinearBridge,
            model -> begin
                x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x1, x2]),
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end,
            model -> begin
                x1, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                x2, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                # x1 * (x2 - 3) <= 0
                MOI.add_constraint(model, x1 * (x2 - 3.0), MOI.LessThan(0.0))
            end;
            cannot_unbridge = true,
        )
    end

    @testset "LessThan with non-zero bound" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.NonlinearBridge,
            model -> begin
                x1, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                x2, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x1, x2]),
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
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
    end
end
