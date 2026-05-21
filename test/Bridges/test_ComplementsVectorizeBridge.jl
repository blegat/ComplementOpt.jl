# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestComplementsVectorizeBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestComplementsVectorizeBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_GreaterThan_Nonnegatives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
        """
        variables: x, y
        x >= 0.0
        y >= 3.0
        [x, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y >= 3.0
        [1.0 * x, 1.0 * y + -3.0] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_LessThan_Nonpositives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
        """
        variables: x, y
        x >= 0.0
        y <= -2.0
        [x, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y <= -2.0
        [1.0 * x, 1.0 * y + 2.0] in $M.ComplementsWithSetType{MOI.Nonpositives}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_EqualTo_Zeros()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ComplementsVectorizeBridge{Float64},
        """
        variables: x, y
        x >= 0.0
        y == 2.0
        [x, y] in $M.ComplementsWithSetType{MOI.EqualTo{Float64}}(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y == 2.0
        [1.0 * x, 1.0 * y + -2.0] in $M.ComplementsWithSetType{MOI.Zeros}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestComplementsVectorizeBridge

TestComplementsVectorizeBridge.runtests()
