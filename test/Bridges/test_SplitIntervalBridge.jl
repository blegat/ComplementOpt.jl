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

end  # TestSplitIntervalBridge

TestSplitIntervalBridge.runtests()
