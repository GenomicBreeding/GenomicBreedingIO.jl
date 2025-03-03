"""
    readdelimited(
        type::Type{Genomes};
        fname::String,
        sep::String = \"\\t\",
        parse_populations_from_entries::Union{Nothing,Function} = nothing,
        verbose::Bool = false,
    )::Genomes

Load a `Genomes` struct from a string-delimited (default=\"\\t\") file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by \"|\"), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.

# Notes:
- Extension name should be '.tsv', '.csv', or '.txt'.
- Header lines and comments are prefixed by '#'.
- There are 2 header lines prefixed by '#', e.g.:
    + header line 1: "chrom,pos,all_alleles,allele,entry_1,entry_2"
    + header line 2: "chrom,pos,all_alleles,allele,population_1,population_1"

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=10, verbose=false);

julia> genomes.entries = [string(genomes.populations[i], "-", genomes.entries[i]) for i in eachindex(genomes.populations)];

julia> fname = writedelimited(genomes);

julia> genomes_reloaded = readdelimited(Genomes, fname=fname);

julia> genomes == genomes_reloaded
true

julia> fname = writedelimited(genomes, include_population_header=false);

julia> genomes_reloaded = readdelimited(Genomes, fname=fname);

julia> unique(genomes_reloaded.populations) == ["Unknown_population"]
true

julia> genomes_reloaded = readdelimited(Genomes, fname=fname, parse_populations_from_entries=x -> split(x, "-")[1]);

julia> genomes == genomes_reloaded
true
```
"""
function readdelimited(
    type::Type{Genomes};
    fname::String,
    sep::String = "\t",
    parse_populations_from_entries::Union{Nothing,Function} = nothing,
    verbose::Bool = false,
)::Genomes
    # type = Genomes; genomes = GBCore.simulategenomes(n=10, sparsity=0.01); fname = writedelimited(genomes); sep = "\t"; verbose = true;
    # parse_populations_from_entries = x -> split(x, "-")[1]
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    # Count the number of lines in the file which are not header lines or comments
    # n_lines::Int64 = countlines(fname) # Only works if there are no comments other than the first 2 header lines
    n_lines::Int64 = 0
    line_counter::Int64 = 0
    file = open(fname, "r")
    for raw_line in eachline(file)
        line_counter += 1
        if (raw_line[1] != '#') && (line_counter > 2) && (length(raw_line) > 0)
            n_lines += 1
        end
    end
    close(file)
    # Read the 2 header lines
    file = open(fname, "r")
    header_1::Vector{String} = split(readline(file), sep)
    header_2::Vector{String} = split(readline(file), sep)
    close(file)
    # Column index of the start of numeric values
    IDX::Int64 = 5
    # Expected header/column names
    expected_colnames = ["chrom", "pos", "all_alleles", "allele"]
    if length(header_1) < IDX
        throw(
            ArgumentError(
                "The genomes file `" *
                fname *
                "` have less than " *
                string(IDX) *
                " columns. We expect in order the following column names: \n\t• " *
                join(expected_colnames, "\n\t• "),
            ),
        )
    end
    idx_mismatch = findall(.!isfuzzymatch.(expected_colnames, header_1[1:(IDX-1)]))
    if length(idx_mismatch) > 0
        throw(
            ArgumentError(
                "The genomes file `" *
                fname *
                "` have the following expected column name mismatches: \n\t• " *
                join(string.(expected_colnames[idx_mismatch], " <=> ", header_1[1:(IDX-1)][idx_mismatch]), "\n\t• "),
            ),
        )
    end
    # If only one header exits, i.e. the header with entriy names and the population names are missing
    with_header_2::Bool = (header_1[1:4] == header_2[1:4])
    header_2 = if with_header_2
        header_2
    else
        population_names = if isnothing(parse_populations_from_entries)
            repeat(["Unknown_population"], length(header_1) - 4)
        else
            try
                parse_populations_from_entries.(header_1[IDX:end])
            catch
                throw(
                    ArgumentError(
                        "Error parsing population names from entry names using the `parse_populations_from_entries` function.",
                    ),
                )
            end
        end
        vcat(header_1[1:4], population_names)
    end
    if (length(header_1) != length(header_2))
        throw(ErrorException("The 2 header lines in the genomes file: '" * fname * "' do not match."))
    end
    # Define the expected dimensions of the Genomes struct
    n::Int64 = length(header_1) - (IDX - 1)
    p::Int64 = if with_header_2
        n_lines
    else
        n_lines + 1
    end
    # Instatiate the output struct
    genomes = Genomes(n = n, p = p)
    genomes.entries = header_1[IDX:end]
    genomes.populations = header_2[IDX:end]
    genomes.mask .= true
    # Check for duplicate entries
    unique_entries::Vector{String} = unique(genomes.entries)
    sorted_entries = sort(genomes.entries)
    duplicated_entries::Vector{String} = []
    for i = 2:length(sorted_entries)
        if sorted_entries[i-1] == sorted_entries[i]
            println(sorted_entries[i-1])
            push!(duplicated_entries, sorted_entries[i-1])
        end
    end
    if length(genomes.entries) > length(unique_entries)
        throw(
            ErrorException(
                string("Duplicate entries in file: '", fname, "' i.e.:\n\t‣ ", join(duplicated_entries, "\n\t‣ ")),
            ),
        )
    end
    # Read the file line by line
    line_counter = 0
    i::Int64 = 0
    file = open(fname, "r")
    allele_frequencies::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(n_lines; desc = "Loading genotype file: ")
    end
    line::Vector{String} = repeat([""], length(header_1))
    for raw_line in eachline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line .= split(raw_line, sep)
        line_counter += 1
        # Skip commented out lines and empty lines, as well as the first 2 lines if we have 2 headers or just the first line if we have only 1 header
        if (
            (line[1][1] != '#') &&
            (length(raw_line) > 0) &&
            ((with_header_2 && (line_counter > 2)) || (!with_header_2 && (line_counter > 1)))
        )
            if length(header_1) != length(line)
                throw(
                    ErrorException(
                        "The header line and line: " *
                        string(line_counter) *
                        " of the genomes file: '" *
                        fname *
                        "' do have the same number of columns.",
                    ),
                )
            end
            i += 1
            chrom = line[1]
            pos = try
                parse(Int64, line[2])
            catch
                throw(
                    ErrorException(
                        "Cannot parse the second column, i.e. position field at line: " *
                        string(line_counter) *
                        " ('" *
                        line[2] *
                        "') of the genomes file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            refalts = line[3]
            allele = line[4]
            genomes.loci_alleles[i] = join([chrom, pos, refalts, allele], "\t")
            # Catch missing allele frequencies and convert to -999 for parsing
            bool_missing =
                (line[IDX:end] .== "missing") .||
                (line[IDX:end] .== "NA") .||
                (line[IDX:end] .== "na") .||
                (line[IDX:end] .== "N/A") .||
                (line[IDX:end] .== "n/a") .||
                (line[IDX:end] .== "")
            idx_missing = findall(bool_missing)
            idx_non_missing = findall(.!bool_missing)
            if length(idx_missing) > 0
                line[collect(IDX:end)[idx_missing]] .= "-999"
            end
            allele_frequencies .= try
                parse.(Float64, line[IDX:end])
            catch
                throw(
                    ErrorException(
                        "Cannot parse columns " *
                        string(IDX) *
                        " to " *
                        string(length(line)) *
                        ", i.e. allele frequencies at line: " *
                        string(line_counter) *
                        " of the genomes file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            # Set -999 allele frequencies to missing
            allele_frequencies[idx_missing] .= missing
            # Catch allele frequency overflows
            if sum((allele_frequencies[idx_non_missing] .< 0.0) .|| (allele_frequencies[idx_non_missing] .> 1.0)) > 0
                throw(
                    OverflowError(
                        "Allele frequencies greater than 1 and/or less than 0 at line: " *
                        string(line_counter) *
                        " of the genomes file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            # Insert the allele frequencies
            genomes.allele_frequencies[:, i] = allele_frequencies
            if verbose
                ProgressMeter.next!(pb)
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    close(file)
    # Checks
    unique_loci_alleles::Vector{String} = unique(genomes.loci_alleles)
    sorted_loci_alleles = sort(genomes.loci_alleles)
    duplicated_loci_alleles::Vector{String} = []
    for i = 2:length(sorted_loci_alleles)
        if sorted_loci_alleles[i-1] == sorted_loci_alleles[i]
            println(sorted_loci_alleles[i-1])
            push!(duplicated_loci_alleles, sorted_loci_alleles[i-1])
        end
    end
    if length(genomes.loci_alleles) > length(unique_loci_alleles)
        throw(
            ErrorException(
                string(
                    "Duplicate loci-allele combinations in file: '",
                    fname,
                    "' at:\n\t‣ ",
                    join(duplicated_loci_alleles, "\n\t‣ "),
                ),
            ),
        )
    end
    if !checkdims(genomes)
        throw(ErrorException("Error loading Genomes struct from the file: '" * fname * "'"))
    end
    # Output
    genomes
end

"""
    writedelimited(
        genomes::Genomes;
        fname::Union{Missing,String} = missing,
        sep::String = \"\\t\",
        include_population_header::Bool = true,
    )::String

Save `Genomes` struct as a string-delimited (default=\"\\t\") file.
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by \"|\"), and the specific allele.
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
function writedelimited(
    genomes::Genomes;
    fname::Union{Missing,String} = missing,
    sep::String = "\t",
    include_population_header::Bool = true,
)::String
    # genomes = Genomes(n=2,p=4); genomes.entries = ["entry_1", "entry_2"]; genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"]; sep::String = "\t"; fname = missing; include_population_header = true;
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
        append!(header_1, genomes.entries)
        header_1[end] *= "\n"
        write(file, join(header_1, sep))
        if include_population_header
            header_2::Vector{String} = ["#chrom", "pos", "all_alleles", "allele"]
            append!(header_2, genomes.populations)
            header_2[end] *= "\n"
            write(file, join(header_2, sep))
        end
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
    readdelimited(type::Type{Phenomes}; fname::String, sep::String = "\\t")::Phenomes

Load a `Phenomes` struct from a string-delimited (default=\"\\t\") file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by \"|\"), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.populations .= "pop1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(10,3); phenomes.phenotypes[1,1] = missing; phenomes.mask .= true;

julia> fname = writedelimited(phenomes);

julia> phenomes_reloaded = readdelimited(Phenomes, fname=fname);

julia> phenomes == phenomes_reloaded
true
```
"""
function readdelimited(type::Type{Phenomes}; fname::String, sep::String = "\t", verbose::Bool = false)::Phenomes
    # type = Phenomes; phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(10,3); phenomes.mask .= true; fname = writedelimited(phenomes); sep = "\t"; verbose = true
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    # Count the number of lines in the file which are not header line or comments
    n_lines::Int64 = 0
    file = open(fname, "r")
    line_counter::Int64 = 0
    for raw_line in eachline(file)
        line_counter += 1
        if (raw_line[1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            n_lines += 1
        end
    end
    close(file)
    # Read the header line
    file = open(fname, "r")
    header::Vector{String} = split(readline(file), sep)
    close(file)
    # Column index of the start of numeric values
    IDX::Int64 = 3
    # Expected header/column names
    expected_colnames = ["entries", "populations"]
    if length(header) < IDX
        throw(
            ArgumentError(
                "The phenomes file `" *
                fname *
                "` have less than " *
                string(IDX) *
                " columns. We expect in order the following column names: \n\t• " *
                join(expected_colnames, "\n\t• "),
            ),
        )
    end
    idx_mismatch = findall(.!isfuzzymatch.(expected_colnames, header[1:(IDX-1)]))
    if length(idx_mismatch) > 0
        throw(
            ArgumentError(
                "The phenomes file `" *
                fname *
                "` have the following expected column name mismatches: \n\t• " *
                join(string.(expected_colnames[idx_mismatch], " <=> ", header[1:(IDX-1)][idx_mismatch]), "\n\t• "),
            ),
        )
    end
    # Define the expected dimensions of the Phenomes struct
    n::Int64 = n_lines
    t::Int64 = length(header) - (IDX - 1)
    # Instatiate the output struct
    phenomes = Phenomes(n = n, t = t)
    phenomes.traits = header[IDX:end]
    phenomes.mask .= true
    # Check for duplicate traits
    unique_traits::Vector{String} = unique(phenomes.traits)
    duplicated_traits::Vector{String} = []
    for trait in unique_traits
        if sum(phenomes.traits .== trait) > 1
            push!(duplicated_traits, trait)
        end
    end
    if length(phenomes.traits) > length(unique_traits)
        throw(
            ErrorException(
                string("Duplicate traits in file: '", fname, "' i.e.:\n\t‣ ", join(duplicated_traits, "\n\t‣ ")),
            ),
        )
    end
    # Read the file line by line
    line_counter = 0
    i::Int64 = 0
    file = open(fname, "r")
    line::Vector{String} = repeat([""], length(header))
    phenotypes::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(n_lines; desc = "Loading phenotype file: ")
    end
    for raw_line in eachline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line .= split(raw_line, sep)
        line_counter += 1
        # Skip commented out lines including the header line and empty lines
        if (line[1][1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            if length(header) != length(line)
                throw(
                    ErrorException(
                        "The header line and line: " *
                        string(line_counter) *
                        " of the genomes file: '" *
                        fname *
                        "' do have the same number of columns.",
                    ),
                )
            end
            i += 1
            phenomes.entries[i] = line[1]
            phenomes.populations[i] = line[2]
            # Catch missing phenotypes and convert to -999 for parsing
            bool_missing =
                (line[IDX:end] .== "missing") .||
                (line[IDX:end] .== "NA") .||
                (line[IDX:end] .== "na") .||
                (line[IDX:end] .== "N/A") .||
                (line[IDX:end] .== "n/a") .||
                (line[IDX:end] .== "")
            idx_missing = findall(bool_missing)
            # idx_non_missing = findall(.!bool_missing)
            if length(idx_missing) > 0
                line[collect(IDX:end)[idx_missing]] .= "-999"
            end
            phenotypes = try
                parse.(Float64, line[IDX:end])
            catch
                throw(
                    ErrorException(
                        "Cannot parse columns " *
                        string(IDX) *
                        " to " *
                        string(length(line)) *
                        ", i.e. numeric trait values at line: " *
                        string(line_counter) *
                        " of the phenomes file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            phenotypes[idx_missing] .= missing
            phenomes.phenotypes[i, :] = phenotypes
            if verbose
                ProgressMeter.next!(pb)
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    close(file)
    # Checks
    unique_entries::Vector{String} = unique(phenomes.entries)
    duplicated_entries::Vector{String} = []
    for locus_allele in unique_entries
        if sum(phenomes.entries .== locus_allele) > 1
            push!(duplicated_entries, locus_allele)
        end
    end
    if length(phenomes.entries) > length(unique_entries)
        throw(
            ErrorException(
                string("Duplicate entries in file: '", fname, "' at:\n\t‣ ", join(duplicated_entries, "\n\t‣ ")),
            ),
        )
    end
    if !checkdims(phenomes)
        throw(ErrorException("Error loading Genomes struct from the file: '" * fname * "'"))
    end
    # Rename for operation symbol_strings into underscores in the trait names
    symbol_strings::Vector{String} = ["+", "-", "*", "/", "%"]
    for i in eachindex(symbol_strings)
        phenomes.traits = replace.(phenomes.traits, symbol_strings[i] => "_")
    end
    # Output
    phenomes
end

"""
    writedelimited(phenomes::Phenomes, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Phenomes` struct as a string-delimited (default=\"\\t\") file. 
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
    readdelimited(type::Type{Trials}; fname::String, sep::String = "\\t")::Trials

Load a `Trials` struct from a string-delimited (default=\"\\t\") file. 
We expect the following 10 identifier columns:
    ‣ years
    ‣ seasons
    ‣ harvests
    ‣ sites
    ‣ entries
    ‣ populations
    ‣ replications
    ‣ blocks
    ‣ rows
    ‣ cols
All other columns are assumed to be numeric phenotype values.

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=10, verbose=false);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, sparsity=0.1, verbose=false);

julia> fname = writedelimited(trials);

julia> trials_reloaded = readdelimited(Trials, fname=fname);

julia> trials == trials_reloaded
true
```
"""
function readdelimited(type::Type{Trials}; fname::String, sep::String = "\t", verbose::Bool = false)::Trials
    # type = Trials; genomes = GBCore.simulategenomes(n=10, verbose=false); trials, _ = GBCore.simulatetrials(genomes=genomes, verbose=false); fname = writedelimited(trials); sep = "\t"; verbose = true;
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    # Count the number of lines in the file which are not header line or comments
    n_lines::Int64 = 0
    file = open(fname, "r")
    line_counter::Int64 = 0
    for raw_line in eachline(file)
        line_counter += 1
        if (raw_line[1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            n_lines += 1
        end
    end
    close(file)
    if verbose
        println(string("Reading a ", n_lines, "-line Trials file."))
    end
    # Read the header line
    file = open(fname, "r")
    header::Vector{String} = split(readline(file), sep)
    close(file)
    # Expected minimum header/column names
    expected_colnames =
        ["years", "seasons", "harvests", "sites", "entries", "populations", "replications", "blocks", "rows", "cols"]
    idx_expected_colnames = []
    for name in expected_colnames
        # name = expected_colnames[4]
        idx = findall(isfuzzymatch.(header, name))
        if length(idx) == 0
            throw(ArgumentError("We cannot find the required identifier column: `" * name * "` in the header line."))
        elseif length(idx) > 1
            levdist = [levenshteindistance(name, header[ix]) for ix in idx]
            idx = idx[findall(levdist .== minimum(levdist))[1]]
        else
            idx = idx[1]
        end
        if (length(idx_expected_colnames) > 0) && (sum(idx_expected_colnames .== idx) > 0)
            throw(
                ArgumentError("The identifier column for: `" * name * "` in the header line may have been misspelled."),
            )
        end
        append!(idx_expected_colnames, idx)
    end
    idx_traits = setdiff(collect(1:length(header)), idx_expected_colnames)
    # Define the expected dimensions of the Trials struct
    n::Int64 = n_lines
    t::Int64 = length(header) - length(expected_colnames)
    # Instatiate the output struct
    trials = Trials(n = n, t = t)
    trials.traits = header[idx_traits]
    # Check for duplicate traits
    unique_traits::Vector{String} = unique(trials.traits)
    duplicated_traits::Vector{String} = []
    for trait in unique_traits
        if sum(trials.traits .== trait) > 1
            push!(duplicated_traits, trait)
        end
    end
    if length(trials.traits) > length(unique_traits)
        throw(
            ErrorException(
                string("Duplicate traits in file: '", fname, "' i.e.:\n\t‣ ", join(duplicated_traits, "\n\t‣ ")),
            ),
        )
    end
    # Read the file line by line
    line_counter = 0
    i = 0
    file = open(fname, "r")
    line::Vector{String} = repeat([""], length(header))
    phenotypes::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(n_lines; desc = "Loading trials file: ")
    end
    for raw_line in eachline(file)
        # raw_line = readline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line .= split(raw_line, sep)
        line_counter += 1
        # Skip commented out lines including the header line and empty lines
        if (line[1][1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            if length(header) != length(line)
                throw(
                    ErrorException(
                        "The header line and line: " *
                        string(line_counter) *
                        " of the genomes file: '" *
                        fname *
                        "' do have the same number of columns.",
                    ),
                )
            end
            i += 1
            trials.years[i] = line[idx_expected_colnames[1]]
            trials.seasons[i] = line[idx_expected_colnames[2]]
            trials.harvests[i] = line[idx_expected_colnames[3]]
            trials.sites[i] = line[idx_expected_colnames[4]]
            trials.entries[i] = line[idx_expected_colnames[5]]
            trials.populations[i] = line[idx_expected_colnames[6]]
            trials.replications[i] = line[idx_expected_colnames[7]]
            trials.blocks[i] = line[idx_expected_colnames[8]]
            trials.rows[i] = line[idx_expected_colnames[9]]
            trials.cols[i] = line[idx_expected_colnames[10]]
            # Catch missing phenotypes and convert to -999 for parsing
            bool_missing =
                (line[idx_traits] .== "missing") .||
                (line[idx_traits] .== "NA") .||
                (line[idx_traits] .== "na") .||
                (line[idx_traits] .== "N/A") .||
                (line[idx_traits] .== "n/a") .||
                (line[idx_traits] .== "")
            idx_missing = findall(bool_missing)
            # Temporarily set missing values to -999 for vectorised string to numeric parsing.
            # Note that the value is irrelevant as we will use the missing values indexes to convert them into missing - so fret not.
            if length(idx_missing) > 0
                line[collect(idx_traits)[idx_missing]] .= "-999"
            end
            phenotypes = try
                parse.(Float64, line[idx_traits])
            catch
                throw(
                    ErrorException(
                        "Cannot parse columns " *
                        string(IDX) *
                        " to " *
                        string(length(line)) *
                        ", i.e. numeric trait values at line: " *
                        string(line_counter) *
                        " of the trials file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            phenotypes[idx_missing] .= missing
            trials.phenotypes[i, :] = phenotypes
            if verbose
                ProgressMeter.next!(pb)
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    close(file)
    if !checkdims(trials)
        throw(ErrorException("Error loading Genomes struct from the file: '" * fname * "'"))
    end
    # Rename for operation symbol_strings into underscores in the trait names
    symbol_strings::Vector{String} = ["+", "-", "*", "/", "%"]
    for i in eachindex(symbol_strings)
        trials.traits = replace.(trials.traits, symbol_strings[i] => "_")
    end
    # Output
    trials
end

"""
    writedelimited(trials::Trials, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Trials` struct as a string-delimited (default=\"\\t\") file. 
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
