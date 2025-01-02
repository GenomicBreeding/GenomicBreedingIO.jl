using GBIO
using Documenter

DocMeta.setdocmeta!(GBIO, :DocTestSetup, :(using GBIO); recursive = true)

makedocs(;
    modules = [GBIO],
    authors = "jeffersonparil@gmail.com",
    sitename = "GBIO.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GBIO.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GBIO.jl", devbranch = "main")
