using GenomicBreedingCore
using GenomicBreedingIO
using Test
using Documenter

Documenter.doctest(GenomicBreedingIO)

@testset "GenomicBreedingIO.jl" begin
    genomes = GenomicBreedingCore.simulategenomes(n = 2, verbose = false)
    fname = writejld2(genomes)
    @test readjld2(Genomes, fname = fname) == genomes
    phenomes = Phenomes(n = 2, t = 2)
    phenomes.entries = ["entry_1", "entry_2"]
    phenomes.traits = ["trait_1", "trait_2"]
    phenomes.mask .= true
    fname = writejld2(phenomes)
    @test readjld2(Phenomes, fname = fname) == phenomes
    trials, _ = simulatetrials(genomes = genomes, verbose = false)
    fname = writejld2(trials)
    @test readjld2(Trials, fname = fname) == trials
    fname = writedelimited(genomes)
    genomes_reloaded = readdelimited(Genomes, fname = fname)
    @test genomes == genomes_reloaded
    fname = writedelimited(phenomes)
    phenomes_reloaded = readdelimited(Phenomes, fname = fname)
    @test phenomes == phenomes_reloaded
    trials, _ = GenomicBreedingCore.simulatetrials(genomes = genomes, verbose = false)
    fname = writedelimited(trials)
    trials_reloaded = readdelimited(Trials, fname = fname)
    @test trials == trials_reloaded
    fname = writevcf(genomes)
    genomes_reloaded = readvcf(fname = fname)
    @test genomes.entries == genomes_reloaded.entries
    @test dimensions(genomes) == dimensions(genomes_reloaded)
end

# Clean-up
for f in readdir()[match.(r"^output-", readdir()) .!= nothing]
    rm(f)
end
