using GBCore
using Test
using Documenter

Documenter.doctest(GBCore)
for f in readdir()[match.(r"^output-", readdir()).!=nothing]
    rm(f)
end

@testset "GBCore.jl" begin
    @test 1 == 1
end
