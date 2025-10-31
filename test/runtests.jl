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

expected_models = Dict(
    Instances.fletcher_leyffer_ex1_model => fletcher_leyffer_ex1_nonlinear_model,
)

function test_model(model_func)
    model = model_func()
    model = Instances.fletcher_leyffer_ex1_model()
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

instances = filter(names(Instances; all = true)) do name
    # The function types start with `#`
    s = String(name)
    endswith(s, "_model") && !startswith(s, "#")
end

@testset "$name" for name in instances
    test_model(getfield(Instances, name))
end

@testset "Relaxation method $(relax)" for relax in [
    ComplementOpt.ScholtesRelaxation(0.0),
    ComplementOpt.FischerBurmeisterRelaxation(1e-8),
    ComplementOpt.LiuFukushimaRelaxation(1e-8),
    ComplementOpt.KanzowSchwarzRelaxation(1e-8),
]
    model = Instances.fletcher_leyffer_ex1_model()
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), relax)

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    # Need to set the bound relaxation explicitly to 0 for LiuFukushimaRelaxation
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.objective_value(model) ≈ 0.5 atol=1e-8
    @test JuMP.value.(model[:z]) ≈ [0.5, 0.5] atol=1e-8
end

