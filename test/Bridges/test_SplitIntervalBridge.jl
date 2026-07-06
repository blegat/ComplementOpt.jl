# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestSplitIntervalBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestSplitIntervalBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_VectorOfVariables_input()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xn <= 0.0
        1.0 * x + -1.0 * xp + -1.0 * xn == 0.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_VectorAffineFunction_input()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        y in Interval(0.0, 1.0)
        [2.0 * x + 1.0, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xn <= 0.0
        2.0 * x + -1.0 * xp + -1.0 * xn == -1.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_bounded_activity()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        x in Interval(-2.0, 4.0)
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        x in Interval(-2.0, 4.0)
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xp <= 4.0
        xn <= 0.0
        xn >= -2.0
        1.0 * x + -1.0 * xp + -1.0 * xn == 0.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_upper_bounded_activity()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        x <= 4.0
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        x <= 4.0
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xp <= 4.0
        xn <= 0.0
        1.0 * x + -1.0 * xp + -1.0 * xn == 0.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_lower_bounded_activity()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        x >= -2.0
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        x >= -2.0
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xn <= 0.0
        xn >= -2.0
        1.0 * x + -1.0 * xp + -1.0 * xn == 0.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_positive_activity()
    # With `x ∈ [1, 4]`, the negative part `xn = min(x, 0)` is fixed to zero.
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        x in Interval(1.0, 4.0)
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        x in Interval(1.0, 4.0)
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xp <= 4.0
        xn <= 0.0
        xn >= 0.0
        1.0 * x + -1.0 * xp + -1.0 * xn == 0.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_bounded_affine_activity()
    # The activity `2x + 1` with `x ∈ [-2, 4]` lies in `[-3, 9]`.
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SplitIntervalBridge,
        """
        variables: x, y
        x in Interval(-2.0, 4.0)
        y in Interval(0.0, 1.0)
        [2.0 * x + 1.0, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """,
        """
        variables: x, y, xp, xn
        x in Interval(-2.0, 4.0)
        y in Interval(0.0, 1.0)
        xp >= 0.0
        xp <= 9.0
        xn <= 0.0
        xn >= -3.0
        2.0 * x + -1.0 * xp + -1.0 * xn == -1.0
        [xp, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        [xn, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_final_touch_bound_modification()
    T = Float64
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        MathOptComplements.Bridges.SplitIntervalBridge{T},
    }(
        inner,
    )
    x = MOI.add_variable(model)
    y, _ = MOI.add_constrained_variable(model, MOI.Interval(zero(T), one(T)))
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x, y]),
        MathOptComplements.ComplementsWithSetType{MOI.Interval{T}}(2),
    )
    F, S = MOI.VariableIndex, MOI.LessThan{T}
    upper_sets(m) =
        [MOI.get(m, MOI.ConstraintSet(), ci) for
        ci in MOI.get(m, MOI.ListOfConstraintIndices{F,S}())]
    # `x` is unbounded so no upper bound is added on `xp`:
    # the only `LessThan` is `xn <= 0`
    MOI.Bridges.final_touch(model)
    @test upper_sets(inner) == [MOI.LessThan(zero(T))]
    # Adding `x <= 4` adds `xp <= 4` at the next final_touch
    c_ub = MOI.add_constraint(model, x, MOI.LessThan(T(4)))
    MOI.Bridges.final_touch(model)
    @test sort(upper_sets(inner); by = s -> s.upper) ==
          [MOI.LessThan(zero(T)), MOI.LessThan(T(4)), MOI.LessThan(T(4))]
    # Tightening to `x <= 2` updates the bound of `xp`
    MOI.set(model, MOI.ConstraintSet(), c_ub, MOI.LessThan(T(2)))
    MOI.Bridges.final_touch(model)
    @test sort(upper_sets(inner); by = s -> s.upper) ==
          [MOI.LessThan(zero(T)), MOI.LessThan(T(2)), MOI.LessThan(T(2))]
    # Removing the bound of `x` removes the bound of `xp`
    MOI.delete(model, c_ub)
    MOI.Bridges.final_touch(model)
    @test upper_sets(inner) == [MOI.LessThan(zero(T))]
    return
end

end  # TestSplitIntervalBridge

TestSplitIntervalBridge.runtests()
