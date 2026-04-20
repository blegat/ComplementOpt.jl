using Test
using JuMP
using ComplementOpt

@testset "FlipSignBridge" begin
    @testset "Nonnegatives → Nonpositives" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    ComplementOpt.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    ComplementOpt.ComplementsWithSetType{MOI.Nonpositives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Nonpositives → Nonnegatives" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    ComplementOpt.ComplementsWithSetType{MOI.Nonpositives}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    ComplementOpt.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "GreaterThan → LessThan" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    ComplementOpt.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    ComplementOpt.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "LessThan → GreaterThan" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.FlipSignBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    ComplementOpt.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.LessThan(-2.0))
                f = MOIU.operate(vcat, Float64, -1.0 * x, 1.0 * y)
                MOI.add_constraint(
                    model,
                    f,
                    ComplementOpt.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end
end
