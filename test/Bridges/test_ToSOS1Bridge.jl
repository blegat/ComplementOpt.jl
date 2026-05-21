# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestToSOS1Bridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestToSOS1Bridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_ToSOS1Bridge()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ToSOS1Bridge,
        """
        variables: x, y
        x >= 0.0
        y >= 0.0
        [x, y] in $M.ComplementsWithSetType{MOI.Nonnegatives}(2)
        """,
        """
        variables: x, y
        x >= 0.0
        y >= 0.0
        [x, y] in SOS1([1.0, 2.0])
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestToSOS1Bridge

TestToSOS1Bridge.runtests()
