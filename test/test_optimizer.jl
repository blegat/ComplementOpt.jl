# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestOptimizer

using JuMP
using Test

import HiGHS
import Ipopt
import MathOptComplements
import MathOptInterface as MOI

include(joinpath(@__DIR__, "instances.jl"))

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_per_constraint_reformulation()
    model = Model()
    @variable(model, x1)
    @variable(model, 0.0 <= y1)
    c1 = @constraint(model, x1 ⟂ y1)
    @variable(model, x2)
    @variable(model, 0.0 <= y2)
    c2 = @constraint(model, x2 ⟂ y2)
    @objective(model, Min, (x1 - 1)^2 + y1^2 + (x2 - 1)^2 + y2^2)
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            Ipopt.Optimizer;
            with_cache_type = Float64,
        );
        with_cache_type = Float64,
    )
    # Default is Scholtes
    MOI.set(
        model,
        MathOptComplements.DefaultComplementarityReformulation(),
        MathOptComplements.ScholtesRelaxation(0.0),
    )
    # Override c1 with FischerBurmeister
    set_attribute(
        c1,
        MathOptComplements.ComplementarityReformulation(),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    @test MOI.supports(
        unsafe_backend(model),
        MathOptComplements.DefaultComplementarityReformulation(),
    )
    b = unsafe_backend(model)
    attr = MathOptComplements.ComplementarityReformulation()
    F = MOI.VectorOfVariables
    S = MOI.Complements
    @test MOI.Bridges.is_bridged(b, S)
    @test MOI.supports_add_constrained_variables(b, S)
    @test !MOI.Bridges.is_variable_bridged(b, S)
    bridge_type = MOI.Bridges.Constraint.concrete_bridge_type(b, F, S)
    @test bridge_type ==
          MathOptComplements.Bridges.SpecifySetTypeBridge{Float64}
    @test MOI.supports(b, attr, bridge_type)
    @test MOI.supports(
        unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test MOI.supports(
        unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test MOI.supports(
        unsafe_backend(model),
        MathOptComplements.ComplementarityReformulation(),
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Complements},
    )
    @test get_attribute(
        c1,
        MathOptComplements.ComplementarityReformulation(),
    ) isa MathOptComplements.FischerBurmeisterRelaxation
    set_attribute(model, "bound_relax_factor", 0.0)
    set_silent(model)
    optimize!(model)
    @test is_solved_and_feasible(model)
    # Test get through the LazyBridgeOptimizer
    lazy = backend(model).optimizer
    @test !MOI.Bridges.is_bridged(lazy, S)
    ci_mapped = first(
        MOI.get(
            lazy,
            MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.Complements}(),
        ),
    )
    @test MOI.get(
        lazy,
        MathOptComplements.ComplementarityReformulation(),
        ci_mapped,
    ) isa MathOptComplements.FischerBurmeisterRelaxation
    @test MOI.get(
        model,
        MathOptComplements.ComplementarityReformulation(),
        c1,
    ) isa MathOptComplements.FischerBurmeisterRelaxation
    @test isnothing(
        MOI.get(model, MathOptComplements.ComplementarityReformulation(), c2),
    )
    return
end

function test_per_constraint_reformulation_after_optimize!()
    model = Model()
    @variable(model, x)
    @variable(model, 0.0 <= y)
    c = @constraint(model, x ⟂ y)
    @objective(model, Min, (x - 1)^2 + y^2)
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            Ipopt.Optimizer;
            with_cache_type = Float64,
        );
        with_cache_type = Float64,
    )
    MOI.set(
        model,
        MathOptComplements.DefaultComplementarityReformulation(),
        MathOptComplements.ScholtesRelaxation(0.0),
    )
    set_attribute(model, "bound_relax_factor", 0.0)
    set_silent(model)
    optimize!(model)
    @test is_solved_and_feasible(model)
    # Scholtes produces a quadratic constraint
    inner = unsafe_backend(unsafe_backend(model).model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarQuadraticFunction{Float64},
            MOI.LessThan{Float64},
        }(),
    ) == 1
    # Change reformulation after first optimize! (bridge.constraints is populated)
    set_attribute(
        c,
        MathOptComplements.ComplementarityReformulation(),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    @test get_attribute(
        c,
        MathOptComplements.ComplementarityReformulation(),
    ) isa MathOptComplements.FischerBurmeisterRelaxation
    optimize!(model)
    @test is_solved_and_feasible(model)
    inner = unsafe_backend(unsafe_backend(model).model)
    # FB produces a nonlinear constraint
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.LessThan{Float64},
        }(),
    ) == 1
    return
end

function test_ComplementsWithSetType_dimension()
    s = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}(4)
    @test MOI.dimension(s) == 4
    return
end

function test_Optimizer_bridge_dispatch()
    # NLP path: VectorAffineFunction in ComplementsWithSetType → NonlinearBridge
    opt_nlp = MathOptComplements.Optimizer(
        () ->
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
    @test MOI.Bridges.is_bridged(opt_nlp, S)
    @test MOI.Bridges.supports_bridging_constrained_variable(opt_nlp, S)
    @test MOI.Bridges.supports_bridging_constraint(
        opt_nlp,
        MOI.VectorAffineFunction{Float64},
        S,
    )
    @test MOI.Bridges.bridge_type(
        opt_nlp,
        MOI.VectorAffineFunction{Float64},
        S,
    ) == MathOptComplements.Bridges.NonlinearBridge{Float64,MOI.Nonnegatives}
    # SOS1 path: VectorOfVariables in ComplementsWithSetType{Zeros} →
    # ToZerosBridge (the activity is complementary to a fixed/free variable, so
    # it is forced to zero).
    opt_sos1 = MathOptComplements.Optimizer(
        HiGHS.Optimizer;
        with_bridge_type = Float64,
    )
    S_zeros = MathOptComplements.ComplementsWithSetType{MOI.Zeros}
    @test MOI.Bridges.bridge_type(opt_sos1, MOI.VectorOfVariables, S_zeros) ==
          MathOptComplements.Bridges.ToZerosBridge{Float64}
    return
end

function test_ComplementarityReformulation_supports_on_bridge_types()
    model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    attr = MathOptComplements.ComplementarityReformulation()
    @test MOI.supports(
        model,
        attr,
        MathOptComplements.Bridges.SplitIntervalBridge,
    )
    @test MOI.supports(model, attr, MathOptComplements.Bridges.FlipSignBridge)
    @test MOI.supports(
        model,
        attr,
        MathOptComplements.Bridges.ComplementsVectorizeBridge,
    )
    @test MOI.supports(model, attr, MathOptComplements.Bridges.NonlinearBridge)
    return
end

function test_ComplementarityReformulation_through_SplitInterval()
    # Stack SplitInterval on top of a model with NonlinearBridge
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    nl_model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    B = MathOptComplements.Bridges.SplitIntervalBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(nl_model)
    x = MOI.add_variable(model)
    y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(2),
    )
    # Default Scholtes → quadratic constraints
    MOI.Bridges.final_touch(model)
    MOI.Bridges.final_touch(nl_model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarQuadraticFunction{Float64},
            MOI.LessThan{Float64},
        }(),
    ) > 0
    # Change to FB and re-reformulate → nonlinear constraints
    attr = MathOptComplements.ComplementarityReformulation()
    MOI.set(
        model,
        attr,
        ci,
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    MOI.Bridges.final_touch(model)
    MOI.Bridges.final_touch(nl_model)
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.LessThan{Float64},
        }(),
    ) > 0
    # Test SplitIntervalBridge metadata
    b_split = MOI.Bridges.bridge(model, ci)
    @test b_split isa MathOptComplements.Bridges.SplitIntervalBridge
    @test length(
        MOI.get(
            b_split,
            MOI.ListOfConstraintIndices{
                MOI.VariableIndex,
                MOI.GreaterThan{Float64},
            }(),
        ),
    ) == 1
    @test length(
        MOI.get(
            b_split,
            MOI.ListOfConstraintIndices{
                MOI.VariableIndex,
                MOI.LessThan{Float64},
            }(),
        ),
    ) == 1
    return
end

function test_ComplementarityReformulation_through_SplitInterval_before_final_touch()
    # Set the per-constraint reformulation *before* `final_touch`, so that
    # `SplitIntervalBridge` stores it and applies it to `lower`/`upper` when it
    # creates them in `final_touch`.
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.NonlinearBridge{Float64}
    nl_model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    B = MathOptComplements.Bridges.SplitIntervalBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(nl_model)
    x = MOI.add_variable(model)
    y, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}(2),
    )
    attr = MathOptComplements.ComplementarityReformulation()
    MOI.set(
        model,
        attr,
        ci,
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
    MOI.Bridges.final_touch(model)
    MOI.Bridges.final_touch(nl_model)
    # Fischer-Burmeister → nonlinear constraints (not the default quadratic).
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarNonlinearFunction,
            MOI.LessThan{Float64},
        }(),
    ) > 0
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{
            MOI.ScalarQuadraticFunction{Float64},
            MOI.LessThan{Float64},
        }(),
    ) == 0
    return
end

function test_ComplementarityReformulation_through_FlipSign()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.FlipSignBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    y, _ = MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.Nonpositives}(2),
    )
    attr = MathOptComplements.ComplementarityReformulation()
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, attr, ci, relax)
    # FlipSign flips Nonpositives → Nonnegatives, check child exists
    S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VectorAffineFunction{Float64},S}(),
    ) == 1
    # Check the attribute was propagated to the child constraint
    b = MOI.Bridges.bridge(model, ci)
    @test MOI.get(inner, attr, b.constraint) isa
          MathOptComplements.FischerBurmeisterRelaxation
    return
end

function test_ComplementarityReformulation_through_ComplementsVectorize()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
    )
    attr = MathOptComplements.ComplementarityReformulation()
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, attr, ci, relax)
    # `ComplementsVectorizeBridge` creates its (shifted) constraint in
    # `final_touch`, applying the reformulation set above.
    MOI.Bridges.final_touch(model)
    # ComplementsVectorize shifts GreaterThan → Nonnegatives, check child exists
    S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
    @test MOI.get(
        inner,
        MOI.NumberOfConstraints{MOI.VectorAffineFunction{Float64},S}(),
    ) == 1
    # Check the attribute was propagated to the child constraint
    b = MOI.Bridges.bridge(model, ci)
    @test MOI.get(inner, attr, b.constraint) isa
          MathOptComplements.FischerBurmeisterRelaxation
    return
end

function test_ComplementarityReformulation_through_ComplementsVectorize_after_final_touch()
    # Same as above but the reformulation is set *after* `final_touch`, so the
    # bridge forwards it to the constraint it already created.
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(3.0))
    ci = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2),
    )
    MOI.Bridges.final_touch(model)
    b = MOI.Bridges.bridge(model, ci)
    @test b.constraint !== nothing
    attr = MathOptComplements.ComplementarityReformulation()
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, attr, ci, relax)
    @test MOI.get(inner, attr, b.constraint) == relax
    return
end

function test_ComplementarityReformulation_through_Vertical_after_final_touch()
    # The left-hand side is an expression, so the constraint goes through
    # `VerticalBridge`, which creates its constraint in `final_touch`. Setting
    # the reformulation afterwards must forward it to that constraint.
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    B = MathOptComplements.Bridges.VerticalBridge{Float64}
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{B}(inner)
    x = MOI.add_variable(model)
    y, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
    ci = MOI.add_constraint(
        model,
        MOI.Utilities.operate(vcat, Float64, 2.0 * x + 3.0, 1.0 * y),
        S(2),
    )
    MOI.Bridges.final_touch(model)
    b = MOI.Bridges.bridge(model, ci)
    @test b.constraint !== nothing
    attr = MathOptComplements.ComplementarityReformulation()
    relax = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    MOI.set(model, attr, ci, relax)
    @test MOI.get(inner, attr, b.constraint) == relax
    return
end

function test_per_constraint_reformulation_with_VerticalBridge()
    # Use an expression LHS so that the constraint goes through VerticalBridge
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    c = @constraint(model, [x + y, y] ∈ MOI.Complements(2))
    @objective(model, Min, x^2 + y^2)
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            Ipopt.Optimizer;
            with_cache_type = Float64,
        );
        with_cache_type = Float64,
    )
    attr = MathOptComplements.ComplementarityReformulation()
    reformulation = MathOptComplements.FischerBurmeisterRelaxation(1e-8)
    set_attribute(c, attr, reformulation)
    set_attribute(model, "bound_relax_factor", 0.0)
    set_silent(model)
    optimize!(model)
    @test MOI.supports(backend(model), attr, typeof(index(c)))
    @test get_attribute(c, attr) == reformulation
    @test is_solved_and_feasible(model)
    return
end

# SpecifySetType → SplitInterval → FlipSign → ToSOS1
function test_bridge_chain()
    # HiGHS does not support ScalarNonlinearFunction, so the Optimizer
    # uses the SOS1 path instead of NonlinearBridge.
    opt = MathOptComplements.Optimizer(
        HiGHS.Optimizer;
        with_bridge_type = Float64,
    )
    # Create a model with an Interval complementarity constraint
    x = MOI.add_variable(opt)
    y, _ = MOI.add_constrained_variable(opt, MOI.Interval(0.0, 1.0))
    ci = MOI.add_constraint(
        opt,
        MOI.VectorOfVariables([x, y]),
        MOI.Complements(2),
    )
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
    return
end

function test_interval_complements_solve_HiGHS()
    # End-to-end solve through the SOS1/MILP path (HiGHS has no native SOS1 nor
    # nonlinear support). The interval bound `0 <= z <= 1e3` makes the second
    # component be shifted, so `VerticalBridge` introduces a slack that must be
    # given finite bounds for `SOS1ToMILPBridge` to apply.
    model = Model(
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(HiGHS.Optimizer; with_bridge_type = Float64),
        ),
    )
    set_silent(model)
    @variable(model, 0 <= z[1:2] <= 1e3)
    @variable(model, z3)
    @objective(model, Min, sum(z) - z3)
    @constraint(model, z3 .<= 4 .* z)
    @constraint(model, z in MOI.Complements(2))
    optimize!(model)
    @test is_solved_and_feasible(model)
    @test isapprox(objective_value(model), 0.0; atol = 1e-6)
    @test isapprox(value.(z), [0.0, 0.0]; atol = 1e-6)
    @test isapprox(value(z3), 0.0; atol = 1e-6)
    return
end

function test_ComplementsWithSetType_added_directly_HiGHS()
    # The user gives the set type directly instead of going through
    # `MOI.Complements`/`SpecifySetTypeBridge`. Because the shift by the lower
    # bound of `x2` is now computed in `final_touch`, the constraint can even be
    # added *before* that bound is set.
    model = Model(
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(HiGHS.Optimizer; with_bridge_type = Float64),
        ),
    )
    set_silent(model)
    @variable(model, x1)
    @variable(model, x2)
    S = MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}
    @constraint(model, [x1, x2] in S(2))
    # Finite bounds are needed by `SOS1ToMILPBridge` (HiGHS has neither native
    # SOS1 nor nonlinear support).
    set_lower_bound(x1, 0.0)
    set_upper_bound(x1, 1e3)
    set_lower_bound(x2, 1.0)
    set_upper_bound(x2, 1e3)
    @constraint(model, x1 + x2 >= 2)
    @objective(model, Min, x1 + 3 * x2)
    optimize!(model)
    @test is_solved_and_feasible(model)
    # The complementarity is `x1 ⟂ (x2 - 1)`, so either `x2 = 1` and `x1 >= 1`
    # (objective `4`), or `x1 = 0` and `x2 >= 2` (objective `6`).
    @test isapprox(objective_value(model), 4.0; atol = 1e-6)
    @test isapprox(value(x1), 1.0; atol = 1e-6)
    @test isapprox(value(x2), 1.0; atol = 1e-6)
    return
end

function test_ComplementsWithSetType_added_directly_Ipopt()
    # Same model as above but the inner solver supports nonlinear constraints,
    # so `NonlinearBridge` reformulates the complementarity instead.
    model = Model(
        () -> MathOptComplements.Optimizer(
            MOI.instantiate(Ipopt.Optimizer; with_cache_type = Float64),
        ),
    )
    set_silent(model)
    @variable(model, x1)
    @variable(model, x2)
    S = MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}
    @constraint(model, [x1, x2] in S(2))
    set_lower_bound(x1, 0.0)
    set_lower_bound(x2, 1.0)
    @constraint(model, x1 + x2 >= 2)
    @objective(model, Min, x1 + 3 * x2)
    optimize!(model)
    @test is_solved_and_feasible(model)
    @test isapprox(objective_value(model), 4.0; atol = 1e-6)
    @test isapprox(value(x1), 1.0; atol = 1e-6)
    @test isapprox(value(x2), 1.0; atol = 1e-6)
    return
end

function test_add_all_bridges_JuMP_GenericModel()
    model = Model(Ipopt.Optimizer)
    MathOptComplements.Bridges.add_all_bridges(model)
    set_silent(model)
    @variable(model, y <= 0)
    @constraint(model, y + 1 ⟂ y)
    optimize!(model)
    @test is_solved_and_feasible(model)
    @test value(y) ≈ -1.0 atol = 1e-7
end

function test_opt_default()
    is_test(f) = startswith("$f", "_test_opt_")
    for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)(MathOptComplements.Optimizer)
    end
    return
end

function test_opt_add_all_bridges()
    is_test(f) = startswith("$f", "_test_opt_")
    function make_opt(inner)
        lazy = MOI.instantiate(inner; with_bridge_type = Float64)
        MathOptComplements.Bridges.add_all_bridges(lazy)
        return lazy
    end
    for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)(make_opt)
    end
    return
end

function _test_opt_lower_bound_nonnegatives(make_opt)
    model = Model()
    @variable(model, x)
    @variable(model, 0.0 <= y)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(MOI.Utilities.Model{Float64}))
    MOI.Utilities.attach_optimizer(model)
    inner = unsafe_backend(model)
    # ComplementsWithSetType is bridged further to nonlinear constraints
    S = MathOptComplements.ComplementsWithSetType{MOI.Nonnegatives}
    @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) ==
          0
    return
end

function _test_opt_lower_bound_greater_than(make_opt)
    model = Model()
    @variable(model, x)
    @variable(model, 3.0 <= y)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(MOI.Utilities.Model{Float64}))
    MOI.Utilities.attach_optimizer(model)
    inner = unsafe_backend(model)
    S = MathOptComplements.ComplementsWithSetType{MOI.GreaterThan{Float64}}
    @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) ==
          0
    return
end

function _test_opt_upper_bound_less_than(make_opt)
    model = Model()
    @variable(model, x)
    @variable(model, y <= 1.0)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(MOI.Utilities.Model{Float64}))
    MOI.Utilities.attach_optimizer(model)
    inner = unsafe_backend(model)
    S = MathOptComplements.ComplementsWithSetType{MOI.LessThan{Float64}}
    @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) ==
          0
    return
end

function _test_opt_range_case_x1_bounded_interval(make_opt)
    model = Model()
    @variable(model, x)
    @variable(model, 0.0 <= y <= 1.0)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(MOI.Utilities.Model{Float64}))
    MOI.Utilities.attach_optimizer(model)
    inner = unsafe_backend(model)
    S = MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}
    @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) ==
          0
    return
end

function _test_opt_range_case_x1_bounded_nonneg(make_opt)
    model = Model()
    @variable(model, 0.0 <= x <= 10.0)
    @variable(model, 0.0 <= y <= 10.0)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(MOI.Utilities.Model{Float64}))
    MOI.Utilities.attach_optimizer(model)
    inner = unsafe_backend(model)
    S = MathOptComplements.ComplementsWithSetType{MOI.Interval{Float64}}
    @test MOI.get(inner, MOI.NumberOfConstraints{MOI.VectorOfVariables,S}()) ==
          0
    return
end

function _test_opt_Vertical(make_opt)
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
    set_optimizer(model, () -> make_opt(Ipopt.Optimizer))
    MOI.Utilities.attach_optimizer(model)
    return
end

function _test_opt_Vertical_errors(make_opt)
    model = Model()
    @variable(model, 0.0 <= x)
    @variable(model, 0.0 <= y)
    @constraint(model, [x, 1.0*y + x] ∈ MOI.Complements(2))
    set_optimizer(model, () -> make_opt(Ipopt.Optimizer))
    @test_throws Exception MOI.Utilities.attach_optimizer(model)
    return
end

_test_opt_Instances(make_opt) = Instances.runtests(make_opt)

function _test_relaxation(relax)
    model = Model()
    @variable(model, 0 <= x <= 1)
    @variable(model, μ)
    @constraint(model, 1 - μ == 0.0)
    @constraint(model, μ ⟂ x)
    vars = [x, μ]
    sol = [0.0, 1.0]
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(
            Ipopt.Optimizer;
            with_cache_type = Float64,
        ),
    )
    set_attribute(
        model,
        MathOptComplements.DefaultComplementarityReformulation(),
        relax,
    )
    set_attribute(model, "bound_relax_factor", 0.0)
    set_silent(model)
    optimize!(model)
    inner = backend(model).optimizer.model.model
    if relax isa MathOptComplements.ScholtesRelaxation
        F = MOI.ScalarQuadraticFunction{Float64}
        G = MOI.ScalarNonlinearFunction
    else
        G = MOI.ScalarQuadraticFunction{Float64}
        F = MOI.ScalarNonlinearFunction
    end
    @test MOI.get(inner, MOI.NumberOfConstraints{F,MOI.LessThan{Float64}}()) > 0
    @test MOI.get(inner, MOI.NumberOfConstraints{G,MOI.LessThan{Float64}}()) ==
          0
    @test is_solved_and_feasible(model)
    @test value.(vars) ≈ sol atol=1e-7
    return
end

function test_relaxation_scholtes()
    return _test_relaxation(MathOptComplements.ScholtesRelaxation(0.0))
end

function test_relaxation_fisher_burmeister()
    return _test_relaxation(
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
    )
end

function test_simple_ncp()
    @testset "$relax" for relax in [
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        @variable(model, y <= 0, start = -1)
        @constraint(model, y + 1 ⟂ y)
        MOI.set(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relax,
        )
        set_attribute(model, "bound_relax_factor", 0.0)
        set_silent(model)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test value(y) ≈ -1.0 atol=1e-7
    end
    return
end

function test_simple_lp_3()
    for relax in [
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model()
        @variable(model, 0.0 <= x[1:2])
        @variable(model, 0.0 <= z[1:2])
        @variable(model, y)
        @constraint(model, -z[1] + y == 0.0)
        @constraint(model, -1.0 - z[2] + y == 0.0)
        @constraint(model, x[1] + x[2] == 1.0)
        @constraint(model, z[1] ⟂ x[1])
        @constraint(model, z[2] ⟂ x[2])
        set_optimizer(
            model,
            () -> MathOptComplements.Optimizer(Ipopt.Optimizer),
        )
        MOI.set(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relax,
        )
        set_attribute(model, "bound_relax_factor", 0.0)
        set_silent(model)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test value([x; z; y]) ≈ [0.0, 1.0, 1.0, 0.0, 1.0] atol=1e-7
    end
    return
end

function test_fletcher_leyffer_ex1_model()
    for relax in [
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        @variable(model, z[1:2])
        set_lower_bound(z[2], 0)
        @objective(model, Min, (z[1] - 1)^2 + z[2]^2)
        @constraint(model, [z[2] - z[1], z[2]] ∈ MOI.Complements(2))
        MOI.set(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relax,
        )
        set_attribute(model, "bound_relax_factor", 0.0)
        set_silent(model)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test objective_value(model) ≈ 0.5 atol=1e-5
        @test value.(model[:z]) ≈ [0.5, 0.5] atol=1e-5
    end
    return
end

function test_reformulation_doesnt_error()
    for relax in [
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        # Case 1: Complementarity defined as lower-bound on RHS
        @variable(model, x1)
        @variable(model, 0.0 <= y1)
        @constraint(model, [x1, y1] ∈ MOI.Complements(2))
        # Case 2: Complementarity defined as upper-bound on RHS
        @variable(model, x2)
        @variable(model, y2 <= 1.0)
        @constraint(model, [x2, y2] ∈ MOI.Complements(2))
        MOI.set(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relax,
        )
        MOI.Utilities.attach_optimizer(model)
    end
    return
end

function test_nonlinear_reformulation()
    model = Model()
    @variable(model, x >= 0.0)
    @variable(model, y >= 0.0)
    @objective(model, Min, x^2 + y^2 - 4*x*y)
    @constraint(model, [sin(x), y] ∈ MOI.Complements(2))
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(MOI.Utilities.Model{Float64}),
    )
    MOI.Utilities.attach_optimizer(model)
    expected = Model()
    @variable(expected, x >= 0.0)
    @variable(expected, y >= 0.0)
    @variable(expected, slack >= 0)
    @objective(expected, Min, x^2 + y^2 - 4*x*y)
    # Build complementarity constraints with nonlinear expression
    @constraint(expected, sin(x) == slack)
    @constraint(expected, slack * y <= 0.0)
    MOI.Bridges._test_structural_identical(
        unsafe_backend(model).model,
        backend(expected),
    )
    return
end

function test_reformulation_fletcher_leyffer_ex1()
    model = Model()
    @variable(model, z[1:2])
    set_lower_bound(z[2], 0)
    @objective(model, Min, (z[1] - 1)^2 + z[2]^2)
    @constraint(model, [z[2] - z[1], z[2]] ∈ MOI.Complements(2))
    set_optimizer(
        model,
        () -> MathOptComplements.Optimizer(MOI.Utilities.Model{Float64}),
    )
    MOI.Utilities.attach_optimizer(model)
    expected = Model()
    @variable(expected, z[1:2])
    @variable(expected, slack >= 0.0)
    set_lower_bound(z[2], 0)
    @objective(expected, Min, (z[1] - 1)^2 + z[2]^2)
    @constraint(expected, z[2] - z[1] == slack)
    @constraint(expected, slack * z[2] <= 0)
    MOI.Bridges._test_structural_identical(
        unsafe_backend(model).model,
        backend(expected),
    )
    return
end

function test_double_sided_bound_reformulations()
    atol = 1e-3
    @testset "$relaxation" for relaxation in Any[
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        # Does not support double sided bounds
        # MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        set_attribute(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relaxation,
        )
        @variable(model, 1 <= x <= 2, start = 1)
        @variable(model, -1 <= y <= 3, start = 3)
        @constraint(model, y - 1 ⟂ x)
        set_silent(model)
        @objective(model, Min, x - 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 1; atol)
        @test isapprox(value(y), 3; atol)
        @objective(model, Min, x + 0.01 * y)
        set_start_value(x, 1)
        set_start_value(y, 1)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 1; atol)
        @test isapprox(value(y), 1; atol)
        @objective(model, Max, x - 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 2; atol)
        @test isapprox(value(y), -1; atol)
        @objective(model, Max, x + 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 2; atol)
        @test isapprox(value(y), 1; atol)
    end
    return
end

function test_greater_than_reformulations()
    atol = 1e-3
    @testset "$relaxation" for relaxation in Any[
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        set_attribute(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relaxation,
        )
        @variable(model, 1 <= x, start = 1)
        @variable(model, -1 <= y <= 3, start = 3)
        @constraint(model, y - 1 ⟂ x)
        set_silent(model)
        @objective(model, Min, x - 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 1; atol)
        @test isapprox(value(y), 3; atol)
        @objective(model, Min, x + 0.01 * y)
        set_start_value(x, 1)
        set_start_value(y, 1)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 1; atol)
        @test isapprox(value(y), 1; atol)
        @objective(model, Max, x)
        optimize!(model)
        @test !is_solved_and_feasible(model)
    end
    return
end

function test_less_than_reformulations()
    atol = 1e-3
    @testset "$relaxation" for relaxation in Any[
        MathOptComplements.ScholtesRelaxation(0.0),
        MathOptComplements.FischerBurmeisterRelaxation(1e-8),
        MathOptComplements.LiuFukushimaRelaxation(1e-8),
        MathOptComplements.KanzowSchwarzRelaxation(1e-8),
    ]
        model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer))
        set_attribute(
            model,
            MathOptComplements.DefaultComplementarityReformulation(),
            relaxation,
        )
        @variable(model, x <= 2, start = 1)
        @variable(model, -1 <= y <= 3, start = 3)
        @constraint(model, y - 1 ⟂ x)
        set_silent(model)
        @objective(model, Max, x - 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 2; atol)
        @test isapprox(value(y), -1; atol)
        @objective(model, Max, x + 0.01 * y)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test isapprox(value(x), 2; atol)
        @test isapprox(value(y), 1; atol)
        @objective(model, Min, x)
        optimize!(model)
        @test !is_solved_and_feasible(model)
    end
    return
end

end  # TestOptimizer

TestOptimizer.runtests()
