using Test
using JuMP
using MathOptComplements

@testset "SplitIntervalBridge" begin
    @testset "VectorOfVariables input" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.SplitIntervalBridge,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([x, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(2),
                )
            end,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
                xp, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                # x == xp + xn  →  x - xp - xn == 0
                MOI.add_constraint(model, 1.0 * x - 1.0 * xp - 1.0 * xn, MOI.EqualTo(0.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([xp, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([xn, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "VectorAffineFunction input" begin
        MOI.Bridges.runtests(
            MathOptComplements.Bridges.SplitIntervalBridge,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
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
                    MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(2),
                )
            end,
            model -> begin
                x = MOI.add_variable(model)
                y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
                xp, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
                # 2x + 1 == xp + xn  →  2x - xp - xn == -1
                MOI.add_constraint(model, 2.0 * x - 1.0 * xp - 1.0 * xn, MOI.EqualTo(-1.0))
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([xp, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
                )
                MOI.add_constraint(
                    model,
                    MOI.VectorOfVariables([xn, y]),
                    MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end
end
