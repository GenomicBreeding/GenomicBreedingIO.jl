using GBCore
using GBIO
using Test
using Documenter

Documenter.doctest(GBIO)

@testset "GBIO.jl" begin
    genomes = GBCore.simulategenomes(n = 2, verbose = false)
    fname = writeJLD2(genomes)
    @test readJLD2(Genomes, fname) == genomes
    phenomes = Phenomes(n = 2, t = 2)
    phenomes.entries = ["entry_1", "entry_2"]
    phenomes.traits = ["trait_1", "trait_2"]
    phenomes.mask .= true
    fname = writeJLD2(phenomes)
    @test readJLD2(Phenomes, fname) == phenomes
    trials, _ = simulatetrials(genomes = genomes, verbose = false)
    fname = writeJLD2(trials)
    @test readJLD2(Trials, fname) == trials
    fname = writedelimited(genomes)
    genomes_reloaded = readdelimited(Genomes, fname = fname)
    @test genomes == genomes_reloaded
    fname = writedelimited(phenomes)
    phenomes_reloaded = readdelimited(Phenomes, fname = fname)
    @test phenomes == phenomes_reloaded
    trials, _ = GBCore.simulatetrials(genomes = genomes, verbose = false)
    fname = writedelimited(trials)
    trials_reloaded = readdelimited(Trials, fname = fname)
    @test trials == trials_reloaded
    fname = writevcf(genomes)
    genomes_reloaded = readvcf(fname = fname)
    @test genomes.entries == genomes_reloaded.entries
    @test dimensions(genomes) == dimensions(genomes_reloaded)
end

# Clean-up
for f in readdir()[match.(r"^output-", readdir()).!=nothing]
    rm(f)
end
