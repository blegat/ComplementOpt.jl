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

# Variable in the left-hand-side should not have two bounds
function test_vertical_mispecified_1()
    model = Model()
    @variable(model, 0.0 <= x <= 1.0)
    @variable(model, 0.0 <= y)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    return model
end

# Variable in the right-hand-side should not be an expression
function test_vertical_mispecified_2()
    model = Model()
    @variable(model, 0.0 <= x)
    @variable(model, 0.0 <= y)
    @constraint(model, [x, 1.0*y + x] ∈ MOI.Complements(2))
    return model
end

function test_vertical_formulation()
    model = Model()
    # Case 1: LHS is already a variable (do nothing)
    @variable(model, x1)
    @variable(model, 0.0 <= y1)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    # Case 2: RHS is unbounded (convert LHS to equality)
    @variable(model, x2)
    @variable(model, y2)
    @constraint(model, [1.0*x2, y2] ∈ MOI.Complements(2))
    # Case 3: LHS is a ScalarAffineFunction with a single variable
    @variable(model, x3)
    @variable(model, 0.0 <= y3)
    @constraint(model, [1.0*x3, y3] ∈ MOI.Complements(2))
    return model
end

function test_nonlinear_reformulation()
    model = Model()
    # Case 1: Complementarity defined as lower-bound on RHS
    @variable(model, x1)
    @variable(model, 0.0 <= y1)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    # Case 2: Complementarity defined as upper-bound on RHS
    @variable(model, x2)
    @variable(model, y2 <= 1.0)
    @constraint(model, [x2, y2] ∈ MOI.Complements(2))
    return model
end

# LHS has a non-trivial lower-bound
function test_nonlinear_mispecified_1()
    model = Model()
    @variable(model, 1.0 <= x1)
    @variable(model, 0.0 <= y1)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    return model
end

# RHS is unbounded
function test_nonlinear_mispecified_2()
    model = Model()
    @variable(model, x1)
    @variable(model, y1)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    return model
end

expected_models =
    Dict(Instances.fletcher_leyffer_ex1_model => fletcher_leyffer_ex1_nonlinear_model)

function test_model(model_func)
    model = model_func()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> ComplementOpt.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)
    if haskey(expected_models, model_func)
        expected_func = expected_models[model_func]
        expected = expected_func()
        MOI.Bridges._test_structural_identical(
            unsafe_backend(model).model,
            backend(expected),
        )
    end
end

function test_nonlinear_expr()
    model = nonlinear_test_model()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> ComplementOpt.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)

    expected = nonlinear_test_reformulated_model()
    MOI.Bridges._test_structural_identical(unsafe_backend(model).model, backend(expected))
end

@testset "Test vertical formulation" begin
    model = test_vertical_formulation()
    set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.Utilities.attach_optimizer(model)

    model = test_vertical_mispecified_2()
    set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    @test_throws Exception MOI.Utilities.attach_optimizer(model)
end

@testset "Nonlinear reformulation $(relax)" for relax in [
    ComplementOpt.ScholtesRelaxation(0.0),
    ComplementOpt.FischerBurmeisterRelaxation(1e-8),
    ComplementOpt.LiuFukushimaRelaxation(1e-8),
    ComplementOpt.KanzowSchwarzRelaxation(1e-8),
]
    model = test_nonlinear_reformulation()
    set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.set(model, ComplementOpt.RelaxationMethod(), relax)
    MOI.Utilities.attach_optimizer(model)

    for test_func in (test_nonlinear_mispecified_1, test_vertical_mispecified_2)
        model = test_func()
        set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
        MOI.set(model, ComplementOpt.RelaxationMethod(), relax)
        @test_throws Exception MOI.Utilities.attach_optimizer(model)
    end
end

instances = filter(names(Instances; all = true)) do name
    # The function types start with `#`
    s = String(name)
    endswith(s, "_model") && !startswith(s, "#")
end

@testset "$name" for name in instances
    test_model(getfield(Instances, name))
end

@testset "Unit-tests" begin
    test_nonlinear_expr()
end

@testset "Relaxation method $(relax)" for relax in [
    ComplementOpt.ScholtesRelaxation(0.0),
    ComplementOpt.FischerBurmeisterRelaxation(1e-8),
    ComplementOpt.LiuFukushimaRelaxation(1e-8),
    ComplementOpt.KanzowSchwarzRelaxation(1e-8),
]
    model = Instances.fletcher_leyffer_ex1_model()

    JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.set(model, ComplementOpt.RelaxationMethod(), relax)
    # Need to set the bound relaxation explicitly to 0 for LiuFukushimaRelaxation
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
    @test JuMP.objective_value(model) ≈ 0.5 atol=1e-7
    @test JuMP.value.(model[:z]) ≈ [0.5, 0.5] atol=1e-7
end
