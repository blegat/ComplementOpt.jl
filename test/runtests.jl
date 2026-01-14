module Instances

include("instances.jl")

end

using Test
using JuMP
using Ipopt

using ComplementOpt

function fletcher_leyffer_ex1_nonlinear_model()
    model = Model()
    @variable(model, z[1:2])
    @variable(model, slack >= 0.0)
    set_lower_bound(z[2], 0)
    @objective(model, Min, (z[1] - 1)^2 + z[2]^2)
    @constraint(model, z[2] - z[1] == slack)
    @constraint(model, slack * z[2] <= 0)
    return model
end

function nonlinear_test_model()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @objective(model, Min, x^2 + y^2 - 4*x*y)
    # Build complementarity constraints with nonlinear expression
    @constraint(model, [sin(x), y] ∈ MOI.Complements(2))
    return model
end

function nonlinear_test_reformulated_model()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @variable(model, slack >= 0)
    @objective(model, Min, x^2 + y^2 - 4*x*y)
    # Build complementarity constraints with nonlinear expression
    @constraint(model, sin(x) == slack)
    @constraint(model, slack * y <= 0.0)
    return model
end

expected_models =
    Dict(Instances.fletcher_leyffer_ex1_model => fletcher_leyffer_ex1_nonlinear_model)

function test_model(model_func)
    model = model_func()
    set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    name = Symbol(model_func)
    @test JuMP.is_solved_and_feasible(model)
    if haskey(Instances.MACMPEC_SOLUTIONS, name)
        @test JuMP.objective_value(model) ≈ Instances.MACMPEC_SOLUTIONS[name] rtol=1e-4 atol=1e-4
    end
end

function test_nonlinear_expr(original_model, reformulated_model)
    model = original_model()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> ComplementOpt.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)
    expected = reformulated_model()
    MOI.Bridges._test_structural_identical(unsafe_backend(model).model, backend(expected))
end

instances = filter(names(Instances; all = true)) do name
    # The function types start with `#`
    s = String(name)
    endswith(s, "_model") && !startswith(s, "#")
end

@testset "$name" for name in instances
    test_model(getfield(Instances, name))
end

@testset "Test reformulation for $original_model" for (
    original_model,
    reformulated_model,
) in [
    (nonlinear_test_model, nonlinear_test_reformulated_model),
    (Instances.fletcher_leyffer_ex1_model, fletcher_leyffer_ex1_nonlinear_model),
]
    test_nonlinear_expr(original_model, reformulated_model)
end

@testset "Relaxation method $(relax)" for relax in [
    ComplementOpt.ScholtesRelaxation(0.0),
    ComplementOpt.FischerBurmeisterRelaxation(1e-8),
    ComplementOpt.LiuFukushimaRelaxation(1e-8),
    ComplementOpt.KanzowSchwarzRelaxation(1e-8),
]
    model = Instances.fletcher_leyffer_ex1_model()

    JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    # Need to set the bound relaxation explicitly to 0 for LiuFukushimaRelaxation
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.objective_value(model) ≈ 0.5 atol=1e-7
    @test JuMP.value.(model[:z]) ≈ [0.5, 0.5] atol=1e-7
end
