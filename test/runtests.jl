module Instances

include("instances.jl")

end

using Test
using JuMP

using ComplementOpt

function fletcher_leyffer_ex1_nonlinear_model()
    model = Model()
    @variable(model, z[1:2])
    @variable(model, slack)
    set_lower_bound(z[2], 0)
    @objective(model, Min, (z[1] - 1)^2 + z[2]^2)
    @constraint(model, z[2] - z[1] >= 0)
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
    @variable(model, slack)
    @objective(model, Min, x^2 + y^2 - 4*x*y)
    # Build complementarity constraints with nonlinear expression
    @constraint(model, sin(x) == slack)
    @constraint(model, slack * y <= 0.0)
    return model
end

expected_models = Dict(
    Instances.fletcher_leyffer_ex1_model => fletcher_leyffer_ex1_nonlinear_model,
)

function test_model(model_func)
    model = model_func()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> ComplementOpt.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)
    MOI.Utilities.attach_optimizer(backend(model).optimizer.model)
    if haskey(expected_models, model_func)
        expected_func = expected_models[model_func]
        expected = expected_func()
        MOI.Bridges._test_structural_identical(
            unsafe_backend(model).optimizer,
            backend(expected),
        )
    end
end

function test_nonlinear_expr()
    model = nonlinear_test_model()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> ComplementOpt.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)
    MOI.Utilities.attach_optimizer(backend(model).optimizer.model)

    expected = nonlinear_test_reformulated_model()
    MOI.Bridges._test_structural_identical(
        unsafe_backend(model).optimizer,
        backend(expected),
    )
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


