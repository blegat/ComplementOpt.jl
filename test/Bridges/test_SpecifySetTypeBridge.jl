# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestSpecifySetTypeBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestSpecifySetTypeBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_lower_bound_Nonnegatives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        x >= 0.0
        y >= 0.0
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y >= 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_lower_bound_GreaterThan()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        x >= 0.0
        y >= 3.0
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y >= 3.0
        [x, y] in $M.ComplementsWithSetType{MOI.GreaterThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_upper_bound_Nonpositives()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        y <= 0.0
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        x <= 0.0
        y <= 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Nonpositives}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_upper_bound_LessThan()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        y <= -2.0
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        x <= 0.0
        y <= -2.0
        [x, y] in $M.ComplementsWithSetType{MOI.LessThan{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_range_Interval()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        y in Interval(0.0, 1.0)
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        y in Interval(0.0, 1.0)
        [x, y] in $M.ComplementsWithSetType{MOI.Interval{Float64}}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_Free_variable_Zeros()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.SpecifySetTypeBridge,
        """
        variables: x, y
        [x, y] in MOI.Complements(2)
        """,
        """
        variables: x, y
        1.0 * x == 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Zeros}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestSpecifySetTypeBridge

TestSpecifySetTypeBridge.runtests()
