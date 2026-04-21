using Test
using JuMP
using ComplementOpt

@testset "ToSOS1Bridge" begin
    MOI.Bridges.runtests(
        ComplementOpt.Bridges.ToSOS1Bridge,
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
            MOI.add_constraint(
                model,
                MOI.VectorOfVariables([x, y]),
                MOI.SOS1([1.0, 2.0]),
            )
        end;
        cannot_unbridge = true,
    )
end
