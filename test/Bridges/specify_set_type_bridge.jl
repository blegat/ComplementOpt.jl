using Test
using JuMP
using MathOptComplements

@testset "SpecifySetTypeBridge" begin
    @testset "Lower bound (Nonnegatives)" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.SpecifySetTypeBridge,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MOI.Complements(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Lower bound (GreaterThan)" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.SpecifySetTypeBridge,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MOI.Complements(2),
                )
            end,
            model -> begin
                x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Upper bound (Nonpositives)" begin
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
    end

    @testset "Upper bound (LessThan)" begin
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
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "Range (Interval)" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.SpecifySetTypeBridge,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MOI.Complements(2),
                )
            end,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end
end
