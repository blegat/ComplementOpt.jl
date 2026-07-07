# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Test

is_test(f) = startswith(f, "test_") && endswith(f, ".jl")

@testset "$root" for (root, dirs, files) in walkdir(@__DIR__)
    @testset "$file" for file in filter(is_test, files)
        include(joinpath(root, file))
    end
end

# Integration test that is not named `test_*.jl` so that it is not picked up by
# the `walkdir` above; it is included explicitly here.
@testset "Bridges/lazy.jl" begin
    include(joinpath(@__DIR__, "Bridges", "lazy.jl"))
end
