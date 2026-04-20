using Test

files_to_exclude = ["runtests.jl"]
for file in readdir(@__DIR__)
    if isdir(joinpath(@__DIR__, file))
        include(joinpath(@__DIR__, file, "runtests.jl"))
    end
    if !endswith(file, ".jl") || any(isequal(file), files_to_exclude)
        continue
    end
    @testset "$(file)" begin
        include(joinpath(@__DIR__, file))
    end
end
