using GenomicBreedingIO
using Documenter

DocMeta.setdocmeta!(GenomicBreedingIO, :DocTestSetup, :(using GenomicBreedingIO); recursive = true)

makedocs(;
    modules = [GenomicBreedingIO],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreedingIO.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreedingIO.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreedingIO.jl", devbranch = "main")
