using Pkg
Pkg.activate(".")
try
    Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
catch
    nothing
end
using GBIO
using GBCore
using FileIO
using JLD2
using GZip
using Dates
using ProgressMeter
