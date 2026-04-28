using Test
using JuMP
using MathOptComplements

@testset "FlipSignBridge" begin
    @testset "Nonnegatives → Nonpositives" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.FlipSignBridge,
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
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Nonpositives → Nonnegatives" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "GreaterThan → LessThan" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "LessThan → GreaterThan" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end
end
