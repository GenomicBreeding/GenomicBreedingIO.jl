"""
    readjld2(type::Type; fname::String)::Type

Load a core (`Genomes`, `Phenomes`, and `Trials`), simulation (`SimulatedEffects`), and model (`TEBV`) struct from a JLD2 file.

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=2, verbose=false);

julia> fname = writejld2(genomes);

julia> readjld2(Genomes, fname=fname) == genomes
true

julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> fname = writejld2(phenomes);

julia> readjld2(Phenomes, fname=fname) == phenomes
true

julia> trials, _ = simulatetrials(genomes=genomes, verbose=false);

julia> fname = writejld2(trials);

julia> readjld2(Trials, fname=fname) == trials
true
```
"""
function readjld2(type::Type{T}; fname::String)::T where {T<:AbstractGB}
    if !isfile(fname)
        throw(ArgumentError("JLD2 file: " * fname * " does not exist."))
    end
    d = load(fname)
    struct_name::String = collect(keys(d))[1]
    x = d[struct_name]
    if !checkdims(x)
        throw(DimensionMismatch(struct_name * " struct from the JLD2 file: " * fname * " is corrupted."))
    end
    return x
end

"""
    writejld2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects}; fname::Union{Missing,String} = missing)::String

Save core (`Genomes`, `Phenomes`, and `Trials`), simulation (`SimulatedEffects`), and model (`TEBV`) structs
as JLD2,  a heirarchical data format version 5 (HDF5) - compatible format.
Note that the extension name should be '.jld2'.

## Examples
```jldoctest; setup = :(using GBCore, GBIO, JLD2)
julia> genomes = GBCore.simulategenomes(n=2, verbose=false);

julia> writejld2(genomes, fname="test_genomes.jld2")
"test_genomes.jld2"

julia> genomes_reloaded = load("test_genomes.jld2");

julia> genomes_reloaded[collect(keys(genomes_reloaded))[1]] == genomes
true

julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> writejld2(phenomes, fname="test_phenomes.jld2")
"test_phenomes.jld2"

julia> phenomes_reloaded = load("test_phenomes.jld2");

julia> phenomes_reloaded[collect(keys(phenomes_reloaded))[1]] == phenomes
true

julia> trials, _ = simulatetrials(genomes=genomes, verbose=false);

julia> writejld2(trials, fname="test_trials.jld2")
"test_trials.jld2"

julia> trials_reloaded = load("test_trials.jld2");

julia> trials_reloaded[collect(keys(trials_reloaded))[1]] == trials
true

julia> simulated_effects = SimulatedEffects();

julia> writejld2(simulated_effects, fname="test_simulated_effects.jld2")
"test_simulated_effects.jld2"

julia> simulated_effects_reloaded = load("test_simulated_effects.jld2");

julia> simulated_effects_reloaded[collect(keys(simulated_effects_reloaded))[1]] == simulated_effects
true

julia> trials, _simulated_effects = GBCore.simulatetrials(genomes = GBCore.simulategenomes(n=10, verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false);

julia> tebv = analyse(trials, max_levels=50, verbose=false);

julia> writejld2(tebv, fname="test_tebv.jld2")
"test_tebv.jld2"

julia> tebv_reloaded = load("test_tebv.jld2");

julia> tebv_reloaded[collect(keys(tebv_reloaded))[1]] == tebv
true
```
"""
function writejld2(A::AbstractGB; fname::Union{Missing,String} = missing)::String
    # Check input arguments
    if !checkdims(A)
        throw(DimensionMismatch(string(typeof(A)) * " input is corrupted."))
    end
    if ismissing(fname)
        fname = string("output-", string(typeof(A)), "-", Dates.format(now(), "yyyymmddHHMMSS"), ".jld2")
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if split(basename(fname), ".")[end] != "jld2"
            throw(ArgumentError("The extension name should be `jld2`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Save as a JLD2 binary file
    save(fname, Dict(string(typeof(A)) => A), compress = true)
    return fname
end
