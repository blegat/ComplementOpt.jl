module Instances

include("instances.jl")

end

using Test
using JuMP
using Ipopt
using HiGHS

using MathOptComplements

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

function simple_ncp()
    model = Model()
    @variable(model, y <= 0)
    @constraint(model, y + 1 ⟂ y)
    return model, [y], [-1.0]
end

# Solve min_x x subject to 0 <= x <= 1
function simple_lp_1()
    model = Model()
    @variable(model, 0 <= x <= 1)
    @variable(model, μ)
    @constraint(model, 1 - μ == 0.0)
    @constraint(model, μ ⟂ x)
    return model, [x, μ], [0.0, 1.0]
end

# Solve min_x -x subject to 0 <= x <= 1
function simple_lp_2()
    model = Model()
    @variable(model, 0 <= x <= 1)
    @variable(model, μ)
    @constraint(model, -1 - μ == 0.0)
    @constraint(model, μ ⟂ x)
    return model, [x, μ], [1.0, -1.0]
end

# Solve min -x2 subject to (x1, x2) >= 0; x1 + x2 = 1
function simple_lp_3()
    model = Model()
    @variable(model, 0.0 <= x[1:2])
    @variable(model, 0.0 <= z[1:2])
    @variable(model, y)
    @constraint(model, -z[1] + y == 0.0)
    @constraint(model, -1.0 - z[2] + y == 0.0)
    @constraint(model, x[1] + x[2] == 1.0)
    @constraint(model, z[1] ⟂ x[1])
    @constraint(model, z[2] ⟂ x[2])
    return model, [x; z; y], [0.0, 1.0, 1.0, 0.0, 1.0]
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

# LHS has a non-trivial upper-bound
function test_nonlinear_mispecified_2()
    model = Model()
    @variable(model, x1 <= 1.0)
    @variable(model, y1 <= 0.0)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    return model
end

# RHS is unbounded
function test_nonlinear_mispecified_3()
    model = Model()
    @variable(model, x1)
    @variable(model, y1)
    @constraint(model, [x1, y1] ∈ MOI.Complements(2))
    return model
end

expected_models =
    Dict(Instances.fletcher_leyffer_ex1_model => fletcher_leyffer_ex1_nonlinear_model)

# The tests that do not set `DefaultComplementarityReformulation` should work
# both with `MathOptComplements.Optimizer` and with a `MOI.Bridges.LazyBridgeOptimizer`
# in which all `MathOptComplements` bridges have been registered via `add_all_bridges`.
const OPTIMIZER_FACTORIES = [
    ("MathOptComplements.Optimizer", inner -> MathOptComplements.Optimizer(inner)),
    (
        "LazyBridgeOptimizer + add_all_bridges",
        inner -> begin
            lazy = MOI.Bridges.full_bridge_optimizer(inner, Float64)
            MathOptComplements.Bridges.add_all_bridges(lazy)
            return lazy
        end,
    ),
]

# Some instances are flaky in CI around the reference optimum, so we loosen
# the tolerance on a per-instance basis.
const MACMPEC_TOLERANCES = Dict(:water_net_model => (rtol = 1e-2, atol = 1e-2))

function test_model(model_func, make_opt)
    model = model_func()
    set_optimizer(model, () -> make_opt(Ipopt.Optimizer()))
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    name = Symbol(model_func)
    @test JuMP.is_solved_and_feasible(model)
    if haskey(Instances.MACMPEC_SOLUTIONS, name)
        tol = get(MACMPEC_TOLERANCES, name, (rtol = 1e-4, atol = 1e-4))
        @test JuMP.objective_value(model) ≈ Instances.MACMPEC_SOLUTIONS[name] rtol=tol.rtol atol=tol.atol
    end
end

function test_nonlinear_expr(original_model, reformulated_model)
    model = original_model()
    inner = MOI.Utilities.Model{Float64}()
    set_optimizer(model, () -> MathOptComplements.Optimizer(inner))
    MOI.Utilities.attach_optimizer(model)
    expected = reformulated_model()
    MOI.Bridges._test_structural_identical(unsafe_backend(model).model, backend(expected))
end

@testset "Test vertical formulation ($(opt_name))" for (opt_name, make_opt) in
                                                       OPTIMIZER_FACTORIES

    model = test_vertical_formulation()
    set_optimizer(model, () -> make_opt(Ipopt.Optimizer()))
    MOI.Utilities.attach_optimizer(model)

    model = test_vertical_mispecified_2()
    set_optimizer(model, () -> make_opt(Ipopt.Optimizer()))
    @test_throws Exception MOI.Utilities.attach_optimizer(model)
end

instances = filter(names(Instances; all = true)) do name
    # The function types start with `#`
    s = String(name)
    endswith(s, "_model") && !startswith(s, "#")
end

@testset "$(opt_name): $name" for (opt_name, make_opt) in OPTIMIZER_FACTORIES,
    name in instances

    test_model(getfield(Instances, name), make_opt)
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

@testset "Relaxation method: $(relax)" for relax in [
    MathOptComplements.ScholtesRelaxation(0.0),
    MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    MathOptComplements.LiuFukushimaRelaxation(1e-8),
    MathOptComplements.KanzowSchwarzRelaxation(1e-8),
]
    @testset "Test reformulation" begin
        model = test_nonlinear_reformulation()
        set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
        MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), relax)
        MOI.Utilities.attach_optimizer(model)

        for test_func in (test_nonlinear_mispecified_1, test_nonlinear_mispecified_2)
            model = test_func()
            set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
            MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), relax)
            @test_throws Exception MOI.Utilities.attach_optimizer(model)
        end
    end

    @testset "Solve Fletcher-Leyffer Ex1 problem" begin
        model = Instances.fletcher_leyffer_ex1_model()
        JuMP.set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
        MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), relax)
        JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
        JuMP.set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.is_solved_and_feasible(model)
        @test JuMP.objective_value(model) ≈ 0.5 atol=1e-7
        @test JuMP.value.(model[:z]) ≈ [0.5, 0.5] atol=1e-7
    end

    @testset "Solve NCP problem $(func)" for func in [simple_ncp, simple_lp_3]
        model, vars, sol = func()
        JuMP.set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
        MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), relax)
        JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
        JuMP.set_silent(model)
        JuMP.optimize!(model)

        @test JuMP.is_solved_and_feasible(model)
        @test JuMP.value.(vars) ≈ sol atol=1e-7
    end
end

# N.B.: at the moment, mixed-complementarity problems are supported only
# with ScholtesRelaxation and with FischerBurmeisterRelaxation
@testset "Mixed-complementarity problem with $(relax)" for relax in [
    MathOptComplements.ScholtesRelaxation(0.0),
    MathOptComplements.FischerBurmeisterRelaxation(1e-8),
]
    @testset "Solve NCP problem $(func)" for func in [simple_lp_1] #, simple_lp_2]
        model, vars, sol = func()
        JuMP.set_optimizer(
            model,
            () -> MathOptComplements.Optimizer(
                MOI.instantiate(Ipopt.Optimizer, with_cache_type = Float64),
            ),
        )
        MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), relax)
        JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
        JuMP.set_silent(model)
        JuMP.optimize!(model)

        inner = backend(model).optimizer.model.model

        if relax isa MathOptComplements.ScholtesRelaxation
            F = MOI.ScalarQuadraticFunction{Float64}
            G = MOI.ScalarNonlinearFunction
        else
            G = MOI.ScalarQuadraticFunction{Float64}
            F = MOI.ScalarNonlinearFunction
        end
        @test MOI.get(inner, MOI.NumberOfConstraints{F,MOI.LessThan{Float64}}()) > 0
        @test MOI.get(inner, MOI.NumberOfConstraints{G,MOI.LessThan{Float64}}()) == 0

        @test JuMP.is_solved_and_feasible(model)
        @test JuMP.value.(vars) ≈ sol atol=1e-7
    end
end

@testset "Per-constraint reformulation" begin
    model = Model()
    @variable(model, x1)
    @variable(model, 0.0 <= y1)
    c1 = @constraint(model, x1 ⟂ y1)
    @variable(model, x2)
    @variable(model, 0.0 <= y2)
    c2 = @constraint(model, x2 ⟂ y2)
    @objective(model, Min, (x1 - 1)^2 + y1^2 + (x2 - 1)^2 + y2^2)
    JuMP.set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(Ipopt.Optimizer, with_cache_type = Float64),
        ),
        with_cache_type = Float64,
    )
    # Default is Scholtes
    MOI.set(
        model,
        MathOptComplements.DefaultComplementarityReformulation(),
        MathOptComplements.ScholtesRelaxation(0.0),
    )
    # Override c1 with FischerBurmeister
    MOI.set(
        model,
        MathOptComplements.ComplementarityReformulation(),
        c1,
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    @test MOI.supports(
        JuMP.unsafe_backend(model),
        MathOptComplements.DefaultComplementarityReformulation(),
    )
    b = JuMP.unsafe_backend(model)
    attr = MathOptComplements.ComplementarityReformulation()
    F = MOI.VectorOfVariables
    S = MOI.Complements
    @test MOI.Bridges.is_bridged(b, S)
    @test MOI.supports_add_constrained_variables(b, S)
    @test !MOI.Bridges.is_variable_bridged(b, S)
    bridge_type = MOI.Bridges.Constraint.concrete_bridge_type(b, F, S)
    @test bridge_type == MathOptComplements.Bridges.SpecifySetTypeBridge{Float64}
    @test MOI.supports(b, attr, bridge_type)
    @test MOI.supports(
        JuMP.unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test MOI.supports(
        JuMP.unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test MOI.supports(
        JuMP.unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test MOI.get(model, MathOptComplements.ComplementarityReformulation(), c1) isa
          MathOptComplements.FischerBurmeisterRelaxation
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    @test JuMP.is_solved_and_feasible(model)
    # Test get through the LazyBridgeOptimizer
    lazy = JuMP.backend(model).optimizer
    @test !MOI.Bridges.is_bridged(lazy, S)
    ci_mapped = first(
        MOI.get(lazy, MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.Complements}()),
    )
    @test MOI.get(lazy, MathOptComplements.ComplementarityReformulation(), ci_mapped) isa
          MathOptComplements.FischerBurmeisterRelaxation
    @test MOI.get(model, MathOptComplements.ComplementarityReformulation(), c1) isa
          MathOptComplements.FischerBurmeisterRelaxation
    @test isnothing(MOI.get(model, MathOptComplements.ComplementarityReformulation(), c2))
end

@testset "Per-constraint reformulation after optimize!" begin
    model = Model()
    @variable(model, x)
    @variable(model, 0.0 <= y)
    c = @constraint(model, x ⟂ y)
    @objective(model, Min, (x - 1)^2 + y^2)
    JuMP.set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(Ipopt.Optimizer, with_cache_type = Float64),
        ),
        with_cache_type = Float64,
    )
    MOI.set(
        model,
        MathOptComplements.DefaultComplementarityReformulation(),
        MathOptComplements.ScholtesRelaxation(0.0),
    )
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    @test JuMP.is_solved_and_feasible(model)
    # Change reformulation after first optimize! (bridge.constraints is populated)
    MOI.set(
        model,
        MathOptComplements.ComplementarityReformulation(),
        c,
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    @test MOI.get(model, MathOptComplements.ComplementarityReformulation(), c) isa
          MathOptComplements.FischerBurmeisterRelaxation
end

@testset "SpecifySetTypeBridge ($(opt_name))" for (opt_name, make_opt) in
                                                  OPTIMIZER_FACTORIES

    @testset "Lower bound (Nonnegatives)" begin
        model = Model()
        @variable(model, x)
        @variable(model, 0.0 <= y)
        @constraint(model, [x, y] ∈ MOI.Complements(2))
        inner = MOI.Utilities.Model{Float64}()
        set_optimizer(model, () -> make_opt(inner))
        MOI.Utilities.attach_optimizer(model)
        # ComplementsWithSetType is bridged further to nonlinear constraints
        S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
        @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) == 0
    end

    @testset "Lower bound (GreaterThan)" begin
        model = Model()
        @variable(model, x)
        @variable(model, 3.0 <= y)
        @constraint(model, [x, y] ∈ MOI.Complements(2))
        inner = MOI.Utilities.Model{Float64}()
        set_optimizer(model, () -> make_opt(inner))
        MOI.Utilities.attach_optimizer(model)
        S = MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}
        @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) == 0
    end

    @testset "Upper bound (LessThan)" begin
        model = Model()
        @variable(model, x)
        @variable(model, y <= 1.0)
        @constraint(model, [x, y] ∈ MOI.Complements(2))
        inner = MOI.Utilities.Model{Float64}()
        set_optimizer(model, () -> make_opt(inner))
        MOI.Utilities.attach_optimizer(model)
        S = MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}
        @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) == 0
    end

    @testset "Range case (x1 unbounded → Interval)" begin
        model = Model()
        @variable(model, x)
        @variable(model, 0.0 <= y <= 1.0)
        @constraint(model, [x, y] ∈ MOI.Complements(2))
        inner = MOI.Utilities.Model{Float64}()
        set_optimizer(model, () -> make_opt(inner))
        MOI.Utilities.attach_optimizer(model)
        S = MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}
        @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) == 0
    end

    @testset "Range case (x1 bounded nonneg)" begin
        model = Model()
        @variable(model, 0.0 <= x <= 10.0)
        @variable(model, 0.0 <= y <= 10.0)
        @constraint(model, [x, y] ∈ MOI.Complements(2))
        inner = MOI.Utilities.Model{Float64}()
        set_optimizer(model, () -> make_opt(inner))
        MOI.Utilities.attach_optimizer(model)
        S = MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}
        @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) == 0
    end
end

@testset "Per-constraint reformulation with VerticalBridge" begin
    # Use an expression LHS so that the constraint goes through VerticalBridge
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    c = @constraint(model, [x + y, y] ∈ MOI.Complements(2))
    @objective(model, Min, x^2 + y^2)
    JuMP.set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(Ipopt.Optimizer, with_cache_type = Float64),
        ),
    )
    attr = MathOptComplements.ComplementarityReformulation()
    reformulation = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, attr, c, reformulation)
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.set_silent(model)
    JuMP.optimize!(model)
    @test MOI.supports(backend(model), attr, typeof(index(c)))
    @test MOI.get(model, attr, c) == reformulation
    @test JuMP.is_solved_and_feasible(model)
end

@testset "Bridge chain: SpecifySetType → SplitInterval → FlipSign → ToSOS1" begin
    # HiGHS does not support ScalarNonlinearFunction, so the Optimizer
    # uses the SOS1 path instead of NonlinearBridge.
    opt = MathOptComplements.Optimizer(
        MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64),
    )

    # Create a model with an Interval complementarity constraint
    x = MOI.add_variable(opt)
    y, _ = MOI.add_constrained_variable(opt, MOI.Interval(0.0, 1.0))
    ci = MOI.add_constraint(opt, MOI.VectorOfVariables([x, y]), MOI.Complements(2))
    MOI.Bridges.final_touch(opt)

    # Step 1: Complements → SpecifySetTypeBridge
    bridge1 = MOI.Bridges.bridge(opt, ci)
    @test bridge1 isa MathOptComplements.Bridges.SpecifySetTypeBridge
    ci_interval = bridge1.constraints[1]
    @test ci_interval isa MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}},
    }

    # Step 2: Interval → SplitIntervalBridge (not NonlinearBridge!)
    bridge2 = MOI.Bridges.bridge(opt, ci_interval)
    @test bridge2 isa MathOptComplements.Bridges.SplitIntervalBridge

    # Step 3: LessThan part → FlipSignBridge
    bridge3 = MOI.Bridges.bridge(opt, bridge2.upper)
    @test bridge3 isa MathOptComplements.Bridges.FlipSignBridge
end

@testset "add_all_bridges(::JuMP.GenericModel)" begin
    model = Model(Ipopt.Optimizer)
    MathOptComplements.Bridges.add_all_bridges(model)
    set_silent(model)
    @variable(model, y <= 0)
    @constraint(model, y + 1 ⟂ y)
    optimize!(model)
    @test is_solved_and_feasible(model)
    @test value(y) ≈ -1.0 atol = 1e-7
end

include(joinpath(@__DIR__, "Bridges", "runtests.jl"))
