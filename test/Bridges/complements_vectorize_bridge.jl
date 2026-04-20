using Test
using JuMP
using ComplementOpt

@testset "ComplementsVectorizeBridge" begin
    @testset "GreaterThan → Nonnegatives" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.ComplementsVectorizeBridge{Float64},
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
                    ComplementOpt.ComplementsWithSetType{MOI.Nonnegatives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end

    @testset "LessThan → Nonpositives" begin
        MOI.Bridges.runtests(
            ComplementOpt.Bridges.ComplementsVectorizeBridge{Float64},
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
                    ComplementOpt.ComplementsWithSetType{MOI.Nonpositives}(2),
                )
            end;
            cannot_unbridge = true,
        )
    end
end
