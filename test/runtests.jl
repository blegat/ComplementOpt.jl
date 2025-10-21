module Instances

include("instances.jl")

end

using Test
using JuMP

function test_model(model::JuMP.Model)
    @test 1 == 1
end

instances = filter(names(Instances; all = true)) do name
    # The function types start with `#`
    s = String(name)
    endswith(s, "_model") && !startswith(s, "#")
end

@testset "$name" for name in instances
    model = getfield(Instances, name)()
    test_model(model)
end
