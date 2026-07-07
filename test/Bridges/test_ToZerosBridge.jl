# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestToZerosBridge

using Test

import MathOptComplements
import MathOptInterface as MOI

const M = "TestToZerosBridge.MathOptComplements"

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_ToZerosBridge()
    MOI.Bridges.runtests(
        MathOptComplements.Bridges.ToZerosBridge,
        """
        variables: x, y
        [x, y] in $M.ComplementsWithSetType{MOI.Zeros}(2)
        """,
        """
        variables: x, y
        [x] in Zeros(1)
        """;
        cannot_unbridge = true,
    )
    return
end

end  # TestToZerosBridge

TestToZerosBridge.runtests()
