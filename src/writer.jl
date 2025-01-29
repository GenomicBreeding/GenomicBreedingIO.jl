"""
    writeJLD2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects}; fname::Union{Missing,String} = missing)::String

Save core (`Genomes`, `Phenomes`, and `Trials`), simulation (`SimulatedEffects`), and model (`TEBV`) structs
as JLD2,  a heirarchical data format version 5 (HDF5) - compatible format.
Note that the extension name should be '.jld2'.

## Examples
```jldoctest; setup = :(using GBCore, GBIO, JLD2)
julia> genomes = GBCore.simulategenomes(n=2, verbose=false);

julia> writeJLD2(genomes, fname="test_genomes.jld2")
"test_genomes.jld2"

julia> genomes_reloaded = load("test_genomes.jld2");

julia> genomes_reloaded[collect(keys(genomes_reloaded))[1]] == genomes
true

julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> writeJLD2(phenomes, fname="test_phenomes.jld2")
"test_phenomes.jld2"

julia> phenomes_reloaded = load("test_phenomes.jld2");

julia> phenomes_reloaded[collect(keys(phenomes_reloaded))[1]] == phenomes
true

julia> trials, _ = simulatetrials(genomes=genomes, verbose=false);

julia> writeJLD2(trials, fname="test_trials.jld2")
"test_trials.jld2"

julia> trials_reloaded = load("test_trials.jld2");

julia> trials_reloaded[collect(keys(trials_reloaded))[1]] == trials
true

julia> simulated_effects = SimulatedEffects();

julia> writeJLD2(simulated_effects, fname="test_simulated_effects.jld2")
"test_simulated_effects.jld2"

julia> simulated_effects_reloaded = load("test_simulated_effects.jld2");

julia> simulated_effects_reloaded[collect(keys(simulated_effects_reloaded))[1]] == simulated_effects
true

julia> trials, _simulated_effects = GBCore.simulatetrials(genomes = GBCore.simulategenomes(n=10, verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false);

julia> tebv = analyse(trials, max_levels=50, verbose=false);

julia> writeJLD2(tebv, fname="test_tebv.jld2")
"test_tebv.jld2"

julia> tebv_reloaded = load("test_tebv.jld2");

julia> tebv_reloaded[collect(keys(tebv_reloaded))[1]] == tebv
true
```
"""
function writeJLD2(A::AbstractGB; fname::Union{Missing,String} = missing)::String
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
    save(fname, Dict(string(typeof(A)) => A))
    return fname
end

"""
    writedelimited(genomes::Genomes, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Genomes` struct as a string-delimited (default=`\t`) file.
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by `|`), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.

## Notes:
- Extension name should be '.tsv', '.csv', or '.txt'.
- Header lines and comments are prefixed by '#'.
- There are 2 header lines prefixed by '#', e.g.:
    + header line 1: "chrom,pos,all_alleles,allele,entry_1,entry_2"
    + header line 2: "chrom,pos,all_alleles,allele,population_1,population_1"

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=2, verbose=false);

julia> writedelimited(genomes, fname="test_genomes.tsv")
"test_genomes.tsv"
```
"""
function writedelimited(genomes::Genomes; fname::Union{Missing,String} = missing, sep::String = "\t")::String
    # genomes = Genomes(n=2,p=4); genomes.entries = ["entry_1", "entry_2"]; genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header lines
        header_1::Vector{String} = ["#chrom", "pos", "all_alleles", "allele"]
        header_2::Vector{String} = ["#chrom", "pos", "all_alleles", "allele"]
        append!(header_1, genomes.entries)
        append!(header_2, genomes.populations)
        header_1[end] *= "\n"
        header_2[end] *= "\n"
        write(file, join(header_1, sep))
        write(file, join(header_2, sep))
        # Rest of the data
        for i = 1:length(genomes.loci_alleles)
            line::Vector{String} = [string(x) for x in split(genomes.loci_alleles[i], "\t")]
            append!(line, [ismissing(x) ? "NA" : string(x) for x in genomes.allele_frequencies[:, i]])
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end

"""
    writedelimited(phenomes::Phenomes, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Phenomes` struct as a string-delimited (default=`\t`) file. 
Each row corresponds to a samples, pools, entries or genotypes.
The first 2 columns correspond to the entry and population names.
The subsequency columns refer to the traits containing the phenotype values of each entry.
Note that the extension name should be '.tsv', '.csv', or '.txt'.

## Notes:
- Extension name should be '.tsv', '.csv', or '.txt'.
- Header line and comments are prefixed by '#'.
- There is 1 header line prefixed by '#', e.g.: "entry,population,trait_1,trait_2,trait_3"

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> writedelimited(phenomes, fname="test_phenomes.tsv")
"test_phenomes.tsv"
```
"""
function writedelimited(phenomes::Phenomes; fname::Union{Missing,String} = missing, sep::String = "\t")::String
    # phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(phenomes)
        throw(DimensionMismatch("Phenomes input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header line
        header::Vector{String} = ["#entries", "populations"]
        append!(header, phenomes.traits)
        header[end] *= "\n"
        write(file, join(header, sep))
        # Rest of the data
        for i in eachindex(phenomes.entries)

            line::Vector{String} = [phenomes.entries[i], phenomes.populations[i]]
            append!(line, [ismissing(x) ? "NA" : string(x) for x in phenomes.phenotypes[i, :]])
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end



"""
    writedelimited(trials::Trials, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Trials` struct as a string-delimited (default=`\t`) file. 
Each row corresponds to a samples, pools, entries or genotypes.
The first 10 columns correspond to:

1. years
2. seasons
3. harvests
4. sites
5. entries
6. populations
7. replications
8. blocks
9. rows
10. cols 

The subsequency columns refer to the traits containing the phenotype values.

## Notes:
- Extension name should be '.tsv', '.csv', or '.txt'.
- Header line and comments are prefixed by '#'.
- There is 1 header line prefixed by '#', e.g.: "years,seasons,harvests, ..., trait_1,tratit_2,trait_3"

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> trials = Trials(n=1, t=2); trials.years = ["year_1"]; trials.seasons = ["season_1"]; trials.harvests = ["harvest_1"]; trials.sites = ["site_1"]; trials.entries = ["entry_1"]; trials.populations = ["population_1"]; trials.replications = ["replication_1"]; trials.blocks = ["block_1"]; trials.rows = ["row_1"]; trials.cols = ["col_1"]; trials.traits = ["trait_1", "trait_2"];

julia> writedelimited(trials, fname="test_trials.tsv")
"test_trials.tsv"
```
"""
function writedelimited(trials::Trials; fname::Union{Missing,String} = missing, sep::String = "\t")::String
    # trials = Trials(n=1, t=2); trials.years = ["year_1"]; trials.seasons = ["season_1"]; trials.harvests = ["harvest_1"]; trials.sites = ["site_1"]; trials.entries = ["entry_1"]; trials.populations = ["population_1"]; trials.replications = ["replication_1"]; trials.blocks = ["block_1"]; trials.rows = ["row_1"]; trials.cols = ["col_1"]; trials.traits = ["trait_1", "trait_2"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(trials)
        throw(DimensionMismatch("Trials input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header line
        header::Vector{String} = [
            "#years",
            "seasons",
            "harvests",
            "sites",
            "entries",
            "populations",
            "replications",
            "blocks",
            "rows",
            "cols",
        ]
        append!(header, trials.traits)
        header[end] *= "\n"
        write(file, join(header, sep))
        # Rest of the data
        for i in eachindex(trials.entries)
            line::Vector{String} = [
                trials.years[i],
                trials.seasons[i],
                trials.harvests[i],
                trials.sites[i],
                trials.entries[i],
                trials.populations[i],
                trials.replications[i],
                trials.blocks[i],
                trials.rows[i],
                trials.cols[i],
            ]
            append!(line, [ismissing(x) ? "NA" : string(x) for x in trials.phenotypes[i, :]])
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end


"""
    writevcf(genomes::Genomes; fname::Union{Missing,String}, ploidy::Int64=0, max_depth::Int64=100, n_decimal_places::Int64=4)::String

Save `Genomes` struct as a variant call format (VCF version 4.2) file.

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes_1 = GBCore.simulategenomes(n=2, verbose=false);

julia> writevcf(genomes_1, fname="test_genomes_1.vcf")
"test_genomes_1.vcf"

julia> genomes_2 = GBCore.simulategenomes(n=2, n_alleles=3, verbose=false);

julia> genomes_2.allele_frequencies = round.(genomes_2.allele_frequencies .* 4) ./ 4;

julia> writevcf(genomes_2, fname="test_genomes_2.vcf", ploidy=4)
"test_genomes_2.vcf"

julia> genomes_3 = GBCore.simulategenomes(n=3, verbose=false);

julia> writevcf(genomes_3, fname="test_genomes_3.vcf", gzip=true)
"test_genomes_3.vcf.gz"
```
"""
function writevcf(
    genomes::Genomes;
    fname::Union{Missing,String} = missing,
    ploidy::Int64 = 0,
    max_depth::Int64 = 100,
    n_decimal_places::Int64 = 4,
    gzip::Bool = false,
)::String
    # genomes = simulategenomes(n_alleles=3, sparsity=0.10); fname = missing; ploidy = 0; max_depth = 100; n_decimal_places = 4; gzip = true;
    # genomes = simulategenomes(n_alleles=3); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 2) ./ 2; fname = missing; ploidy = 2; max_depth = 100; n_decimal_places = 4; gzip = true;
    # genomes = simulategenomes(n_alleles=3); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 4) ./ 4; fname = missing; ploidy = 4; max_depth = 100; n_decimal_places = 4; gzip = true;
    # Check input arguments
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted."))
    end
    if ismissing(fname)
        fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".vcf")
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if split(basename(fname), ".")[end] != "vcf"
            throw(ArgumentError("The extension name should be either `vcf`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Extract locus-allele coordinates
    _chromosomes, _positions, loci_ini_idx, loci_fin_idx = loci(genomes)
    n_alleles = maximum(loci_fin_idx - loci_ini_idx .+ 2)
    # Define the header lines
    header_lines = [
        "##fileformat=VCFv4.2",
        "##fileDate=" * Dates.format(now(), "yyyymmdd"),
        "##source=GBIO.jl-v1.0.0-DEV",
        "##reference=unknown",
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=AD,Number=" * string(n_alleles) * ",Type=Float,Description=\"Allele Depth\">",
        "##FORMAT=<ID=AF,Number=" * string(n_alleles) * ",Type=Float,Description=\"Allele Frequency\">",
        join(vcat(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], genomes.entries), "\t"),
        "",
    ]
    # Write into a new text file
    file = if gzip
        fname = if split(fname, ".")[end] == "vcf"
            replace(fname, ".vcf" => ".vcf.gz")
        end
        GZip.open(fname, "w")
    else
        FileIO.open(fname, "w")
    end
    # Header lines
    if gzip
        GZip.write(file, join(header_lines, "\n"))
    else
        write(file, join(header_lines, "\n"))
    end
    # Rest of the data
    locus_allele_id::Vector{String} = repeat([""], 4)
    line::Vector{String} = repeat([""], length(genomes.entries) + 9)
    for i in eachindex(loci_ini_idx)
        # i = 2
        ini = loci_ini_idx[i]
        fin = loci_fin_idx[i]
        locus_allele_id .= [string(x) for x in split(genomes.loci_alleles[ini], "\t")]
        line[1] = locus_allele_id[1]
        line[2] = locus_allele_id[2]
        line[3] = replace(genomes.loci_alleles[ini], "\t" => "_")
        line[4] = split(locus_allele_id[3], "|")[1]
        line[5] = join(split(locus_allele_id[3], "|")[2:end], ",")
        line[6] = "."
        line[7] = "."
        line[8] = string("NS=", sum(.!ismissing.(genomes.allele_frequencies[:, i])))
        line[9] = "GT:AD:AF"
        # Identify order of the alleles in the vcf output and in the genomes input
        ref = line[4]
        alt = split(line[5], ",")
        alleles = [split(x, "\t")[end] for x in genomes.loci_alleles[ini:fin]]
        append!(alleles, symdiff(split(locus_allele_id[3], "|"), alleles))
        # Extract the allele frequencies including all alleles
        allele_frequencies = genomes.allele_frequencies[:, ini:fin]
        allele_frequencies =
            round.(hcat(allele_frequencies, 1.00 .- sum(allele_frequencies, dims = 2)), digits = n_decimal_places)
        # Set missing data into zero allele frequencies
        allele_frequencies[ismissing.(allele_frequencies)] .= 0.0
        # Make sure the allele frequencies sum up to one per locus
        allele_frequencies = allele_frequencies ./ sum(allele_frequencies, dims = 2)
        # Set NaNs to zeroes. These NaNs resulted from division by zero which in turn is due to missing data.
        allele_frequencies[isnan.(allele_frequencies)] .= 0.0
        # Extract the allele depths
        allele_depths = Int64.(round.(allele_frequencies .* max_depth))
        # Sort according to the order of ref and alts
        idx_col_sort_af_ad = [findall(alleles .== a)[1] for a in vcat(ref, alt)]
        allele_frequencies = allele_frequencies[:, idx_col_sort_af_ad]
        allele_depths = allele_depths[:, idx_col_sort_af_ad]
        # Extract the genotypes
        genotypes = Int64.(allele_frequencies .* ploidy)
        # Define the GT, AD, and AF fields
        gt_tmp = stack(
            [join.(split.(repeat.([alleles[j]], genotypes[:, j]), ""), "/") for j in eachindex(alleles)],
            dims = 2,
        )
        gt_ad_af = repeat([""], length(genomes.entries))
        for j in eachindex(gt_ad_af)
            g = gt_tmp[j, :]
            g = g[g.!=""]
            g = replace.(g, ref => "0")
            for k in eachindex(alt)
                g = replace.(g, alt[k] => string(k))
            end
            if length(g) > 0
                # Sort GT labels so that the values are increasing and then join
                g = split(join(g, "/"), "/")
                sort!(g)
                gt_ad_af[j] = join(g, "/")
            else
                # Missing GT field
                gt_ad_af[j] = "./."
            end
            gt_ad_af[j] = string(gt_ad_af[j], ":", join(allele_depths[j, :], ","))
            gt_ad_af[j] = string(gt_ad_af[j], ":", join(allele_frequencies[j, :], ","))
        end
        line[10:end] = gt_ad_af
        line[end] *= "\n"
        if gzip
            GZip.write(file, join(line, "\t"))
        else
            write(file, join(line, "\t"))
        end
    end
    close(file)
    return fname
end
