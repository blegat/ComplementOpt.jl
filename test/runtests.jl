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
