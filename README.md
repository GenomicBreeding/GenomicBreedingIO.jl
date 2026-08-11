# GenomicBreedingIO

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://genomicbreeding.github.io/GenomicBreedingIO.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://genomicbreeding.github.io/GenomicBreedingIO.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreedingIO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GenomicBreeding/GenomicBreedingIO.jl/actions/workflows/CI.yml?query=branch%3Amain)

Input and output library for GenomicBreeding.jl

## Dev stuff:

### REPL prelude

```shell
julia --project=. -e 'using Pkg; Pkg.instantiate()' # For a fresh Julia installation
julia --project=. --threads=2,1 --load test/interactive_prelude.jl
```

### Format and test

```shell
time julia --project=. --threads=2 -e "using Pkg; Pkg.update()"
time julia --project=. --threads=2  test/cli_tester.jl
```

### Force stable docs on new release:

First create a new release and if the stable don't get updated then try the following:

```shell
TAG=v0.3.0 # should match the release tag you created on Github
git switch main
git push origin :refs/tags/$TAG
git tag -d $TAG
git tag $TAG
git push origin $TAG
git show ${TAG}:Project.toml
git branch --contains $TAG
git commit --allow-empty -m "Rebuild docs"
git push
git switch gh-pages
git pull
git branch --contains v0.3.0
```