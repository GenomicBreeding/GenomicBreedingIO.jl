using Pkg
Pkg.activate(".")
try
    Pkg.add(url = "https://github.com/GenomicBreeding/GenomicBreedingCore.jl")
catch
    nothing
end
using GenomicBreedingIO
using GenomicBreedingCore
using FileIO
using JLD2, CodecZlib
using GZip
using Dates
using ProgressMeter
