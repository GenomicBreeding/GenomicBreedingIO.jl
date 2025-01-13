using Pkg
Pkg.activate(".")
Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
using GBIO
using GBCore
using FileIO
using JLD2
using Dates
using ProgressMeter
