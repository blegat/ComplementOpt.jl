# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestVerticalBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestVerticalBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_Unbounded_RHS_with_constant_in_LHS_Complements()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.VerticalBridge,
        """
        variables: w, x, y, z
        w >= 0.0
        [2.0 * x + 3.0, 1.0 * z, 1.0 * y, 1.0 * w] in MOI.Complements(4)
        """,
        """
        variables: w, x, y, z
        w >= 0.0
        2.0 * x == -3.0
        [z, w] in MOI.Complements(2)
        """;
        cannot_unbridge = true,
    )
    return
end

function test_Expression_LHS_with_constant_ComplementsWithSetType()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.VerticalBridge,
        """
        variables: x, y
        y >= 0.0
        [2.0 * x + 3.0, 1.0 * y] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """,
        """
        variables: x, y, z
        y >= 0.0
        2.0 * x + -1.0 * z == -3.0
        [z, y] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestVerticalBridge

TestVerticalBridge.runtests()
