"""
    readJLD2(type::Type, fname::String = missing)::Type

Load a core (`Genomes`, `Phenomes`, and `Trials`), simulation (`SimulatedEffects`), and model (`TEBV`) struct from a JLD2 file.

## Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=2, sparsity=0.01, verbose=false);

julia> fname = writeJLD2(genomes);

julia> readJLD2(Genomes, fname) == genomes
true

julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> fname = writeJLD2(phenomes);

julia> readJLD2(Phenomes, fname) == phenomes
true

julia> trials = Trials(n=1, t=2); trials.entries = ["entry_1"];

julia> fname = writeJLD2(trials);

julia> readJLD2(Trials, fname) == trials
true
```
"""
function readJLD2(type::Type{T}, fname::String = missing)::T where {T<:AbstractGB}
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
    readdelimited(type::Type{Genomes}; fname::String, sep::String = "\\t")::Genomes

Load a `Genomes` struct from a string-delimited (default=`\\t`) file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by `|`), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=10, sparsity=0.1, verbose=false);

julia> fname = writedelimited(genomes);

julia> genomes_reloaded = readdelimited(Genomes, fname=fname);

julia> genomes == genomes_reloaded
true
```
"""
function readdelimited(type::Type{Genomes}; fname::String, sep::String = "\t", verbose::Bool = false)::Genomes
    # type = Genomes; genomes = GBCore.simulategenomes(n=10, sparsity=0.01); sep = "\t"; fname = writedelimited(genomes); verbose = true
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
    if (length(header_1) != length(header_2)) || (header_1[1:4] != header_2[1:4])
        throw(ErrorException("The 2 header lines in the genomes file: '" * fname * "' do not match."))
    end
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
    # Define the expected dimensions of the Genomes struct
    n::Int64 = length(header_1) - (IDX - 1)
    p::Int64 = n_lines
    # Instatiate the output struct
    genomes = Genomes(n = n, p = p)
    genomes.entries = header_1[IDX:end]
    genomes.populations = header_2[IDX:end]
    genomes.mask .= true
    # Check for duplicate entries
    unique_entries::Vector{String} = unique(genomes.entries)
    duplicated_entries::Vector{String} = []
    for entry in unique_entries
        if sum(genomes.entries .== entry) > 1
            push!(duplicated_entries, entry)
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
    for raw_line in eachline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line = split(raw_line, sep)
        line_counter += 1
        # Skip commented out lines including the first 2 header and empty lines
        if (line[1][1] != '#') && (line_counter > 2) && (length(raw_line) > 0)
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
    duplicated_loci_alleles::Vector{String} = []
    for locus_allele in unique_loci_alleles
        if sum(genomes.loci_alleles .== locus_allele) > 1
            push!(duplicated_loci_alleles, locus_allele)
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
    readdelimited(type::Type{Phenomes}; fname::String, sep::String = "\\t")::Phenomes

Load a `Phenomes` struct from a string-delimited (default=`\\t`) file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by `|`), and the specific allele.
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
    phenotypes::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(n_lines; desc = "Loading phenotype file: ")
    end
    for raw_line in eachline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line = split(raw_line, sep)
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
    readdelimited(type::Type{Trials}; fname::String, sep::String = "\\t")::Trials

Load a `Trials` struct from a string-delimited (default=`\\t`) file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by `|`), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.

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
    # Read the header line
    file = open(fname, "r")
    header::Vector{String} = split(readline(file), sep)
    close(file)
    # Column index of the start of numeric values
    IDX::Int64 = 11
    # Expected header/column names
    expected_colnames =
        ["years", "seasons", "harvests", "sites", "entries", "populations", "replications", "blocks", "rows", "cols"]
    if length(header) < IDX
        throw(
            ArgumentError(
                "The trials file `" *
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
                "The trials file `" *
                fname *
                "` have the following expected column name mismatches: \n\t• " *
                join(string.(expected_colnames[idx_mismatch], " <=> ", header[1:(IDX-1)][idx_mismatch]), "\n\t• "),
            ),
        )
    end
    # Define the expected dimensions of the Trials struct
    n::Int64 = n_lines
    t::Int64 = length(header) - (IDX - 1)
    # Instatiate the output struct
    trials = Trials(n = n, t = t)
    trials.traits = header[IDX:end]
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
    i::Int64 = 0
    file = open(fname, "r")
    phenotypes::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(n_lines; desc = "Loading trials file: ")
    end
    for raw_line in eachline(file)
        # raw_line = readline(file)
        # println(string("i=", i, "; line_counter=", line_counter))
        line = split(raw_line, sep)
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
            trials.years[i] = line[1]
            trials.seasons[i] = line[2]
            trials.harvests[i] = line[3]
            trials.sites[i] = line[4]
            trials.entries[i] = line[5]
            trials.populations[i] = line[6]
            trials.replications[i] = line[7]
            trials.blocks[i] = line[8]
            trials.rows[i] = line[9]
            trials.cols[i] = line[10]
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
    readvcf(;fname::String, field::String = "any", verbose::Bool = false)::Genomes

Load Genomes struct from vcf file

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=10, sparsity=0.1, verbose=false);

julia> fname = writevcf(genomes);

julia> genomes_reloaded = readvcf(fname=fname);

julia> genomes.entries == genomes_reloaded.entries
true

julia> dimensions(genomes) == dimensions(genomes_reloaded)
true

julia> ismissing.(genomes.allele_frequencies) == ismissing.(genomes_reloaded.allele_frequencies)
true
```
"""
function readvcf(; fname::String, field::String = "any", verbose::Bool = false)::Genomes
    # genomes = GBCore.simulategenomes(n=10, sparsity=0.1); fname = writevcf(genomes); field = "any"; verbose = true;
    # genomes = simulategenomes(n_alleles=3, sparsity=0.1); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 4) ./ 4; fname = writevcf(genomes, ploidy=4); field = "GT"; verbose = true;
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    # Count the number of lines in the file which are not header lines or comments
    n_lines::Int64 = 0
    line_counter::Int64 = 0
    file = open(fname, "r")
    for raw_line in eachline(file)
        line_counter += 1
        if (raw_line[1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            n_lines += 1
        end
    end
    close(file)
    # Capture the header containing the names of the entries
    file = open(fname, "r")
    header_line = ""
    for raw_line in eachline(file)
        if !isnothing(match(r"^#CHR", raw_line))
            header_line = raw_line
            break
        end
    end
    close(file)
    header::Vector{String} = split(header_line, "\t")
    # Column index of the start of genotype data
    IDX::Int64 = 10
    # Expected header/column names
    expected_colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    if length(header) < IDX
        throw(
            ArgumentError(
                "The vcf file `" *
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
                "The vcf file `" *
                fname *
                "` have the following expected column name mismatches: \n\t• " *
                join(string.(expected_colnames[idx_mismatch], " <=> ", header[1:(IDX-1)][idx_mismatch]), "\n\t• "),
            ),
        )
    end
    # Extract the names of the entries or samples
    entries = header[IDX:end]
    # Capture the header lines containing the format fields
    file = open(fname, "r")
    format_lines::Vector{String} = []
    for raw_line in eachline(file)
        if !isnothing(match(r"^##FORMAT", raw_line))
            push!(format_lines, raw_line)
        end
        if isnothing(match(r"^#", raw_line))
            break
        end
    end
    close(file)
    # Choose field where priority order starts with AF with the highest priority followed by AD, and finally GT
    if field == "any"
        idx = findall(.!isnothing.(match.(r"ID=AF", format_lines)))
        if length(idx) == 0
            idx = findall(.!isnothing.(match.(r"ID=AD", format_lines)))
        end
        if length(idx) == 0
            idx = findall(.!isnothing.(match.(r"ID=GT", format_lines)))
        end
        if length(idx) == 0
            throw(
                ArgumentError("The input vcf file: `" * fname * "` does not have `AF`, `AD`, or `GT` genotype fields."),
            )
        end
    else
        idx = findall(.!isnothing.(match.(Regex(field), format_lines)))
        if length(idx) == 0
            throw(ArgumentError("The input vcf file: `" * fname * "` does not have `" * field * "` field."))
        end
    end
    # Define the field and the maximum number of alleles per locus
    format_details = split(format_lines[idx[1]], ",")
    field = split(format_details[.!isnothing.(match.(r"ID=", format_details))][1], "=")[end]
    n_alleles = 0
    ploidy = 0 # only for "GT" field - not used for "AD" and "AF" fields
    if field == "GT"
        # Computationally expensive allele counting for field=="GT":
        file = open(fname, "r")
        for raw_line in eachline(file)
            line = split(raw_line, "\t")
            if line[1][1] == '#'
                continue
            end
            idx_field = findall(split(line[9], ":") .== field)
            ψ = length(vcat(split.(split(split(line[IDX], ":")[idx_field[1]], "/"), "|")...))
            a = length(split(line[5], ",")) + 1
            if n_alleles < a
                n_alleles = a
            end
            if ploidy < ψ
                ploidy = ψ
            end
        end
        close(file)
    else
        # Easy n_alleles extraction for "AF" and "AD" fields
        n_alleles = try
            parse(Int64, split(format_details[.!isnothing.(match.(r"Number=", format_details))][1], "=")[end])
        catch
            throw(
                ErrorException(
                    "Error parsing the number of alleles in the `" *
                    field *
                    "` header line of the vcf file: `" *
                    fname *
                    "`",
                ),
            )
        end
    end
    # Define the expected dimensions of the Genomes struct
    n::Int64 = length(entries)
    p::Int64 = n_lines * (n_alleles - 1)
    # Instatiate the output struct
    genomes::Genomes = Genomes(n = n, p = p)
    genomes.entries = entries
    genomes.populations .= "unknown"
    genomes.mask .= true
    # Check for duplicate entries
    unique_entries::Vector{String} = unique(genomes.entries)
    duplicated_entries::Vector{String} = []
    for entry in unique_entries
        if sum(genomes.entries .== entry) > 1
            push!(duplicated_entries, entry)
        end
    end
    if length(genomes.entries) > length(unique_entries)
        throw(
            ErrorException(
                string("Duplicate entries in vcf file: '", fname, "' i.e.:\n\t‣ ", join(duplicated_entries, "\n\t‣ ")),
            ),
        )
    end
    # Read the file line by line
    line_counter = 0
    i::Int64 = 0
    file = open(fname, "r")
    allele_frequencies::Vector{Union{Missing,Float64}} = fill(missing, n)
    if verbose
        pb = ProgressMeter.Progress(p; desc = "Loading genotypes from vcf file: ")
    end
    for raw_line in eachline(file)
        # raw_line = readline(file)
        line = split(raw_line, "\t")
        line_counter += 1
        # Skip commented out lines including the first 2 header
        if line[1][1] != '#'
            if (length(entries) + 9) != length(line)
                throw(
                    ErrorException(
                        "The header line and line: " *
                        string(line_counter) *
                        " of the vcf file: '" *
                        fname *
                        "' do have the same number of columns.",
                    ),
                )
            end
            # Find the field from which we will extract the allele frequencies from
            idx_field = findall(split(line[9], ":") .== field)
            # Skip the line if the field is absent
            if length(idx_field) == 0
                continue
                if verbose
                    println(
                        "Field `" *
                        field *
                        "` absent in line " *
                        string(line_counter) *
                        " of vcf file: `" *
                        fname *
                        "`.",
                    )
                end
            end
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
                        "') of the vcf file: '" *
                        fname *
                        "'.",
                    ),
                )
            end
            ref = line[4]
            alt = split(line[5], ",")
            refalts = vcat([ref], alt)
            if field == "AF"
                afreqs = try
                    parse.(Float64, stack([split(split(x, ":")[idx_field[1]], ",") for x in line[IDX:end]], dims = 1))
                catch
                    throw(
                        ErrorException(
                            "Cannot parse the `" *
                            field *
                            "` field (index=" *
                            string(idx_field[1]) *
                            ") at line " *
                            string(line_counter) *
                            " of the vcf file: `" *
                            fname *
                            "`.",
                        ),
                    )
                end
                idx_missing = findall(sum(afreqs, dims = 2)[:, 1] .== 0.0)
            elseif field == "AD"
                depths = try
                    parse.(Float64, stack([split(split(x, ":")[idx_field[1]], ",") for x in line[IDX:end]], dims = 1))
                catch
                    throw(
                        ErrorException(
                            "Cannot parse the `" *
                            field *
                            "` field (index=" *
                            string(idx_field[1]) *
                            ") at line " *
                            string(line_counter) *
                            " of the vcf file: `" *
                            fname *
                            "`.",
                        ),
                    )
                end
                idx_missing = findall(sum(depths, dims = 2)[:, 1] .== 0.0)
            elseif field == "GT"
                genotype_calls = fill(0, length(line) - 9, ploidy)
                genotype_calls_tmp::Vector{String} = repeat([""], ploidy)
                for (k, x) in enumerate(line[IDX:end])
                    # k, x = 1, line[IDX]
                    genotype_calls_tmp = vcat(split.(split(split(x, ":")[idx_field[1]], "/"), "|")...)
                    # Convert missing GT (i.e. '.') into zeroes
                    replace!(genotype_calls_tmp, "." => "0")
                    genotype_calls[k, :] = parse.(Int64, genotype_calls_tmp)
                end
                idx_missing = findall(sum(genotype_calls, dims = 2)[:, 1] .== 0.0)
            end
            for j in eachindex(refalts)
                # j = 2
                # Skip the reference allele
                if j == 1
                    continue
                end
                allele = refalts[j]
                i += 1
                genomes.loci_alleles[i] = join([chrom, pos, join(refalts, "|"), allele], "\t")
                if field == "AF"
                    # Extract from AF (allele frequencies) field
                    allele_frequencies = afreqs[:, j]
                elseif field == "AD"
                    # Extract from AD (allele depths) field
                    allele_frequencies = depths[:, j] ./ sum(depths, dims = 2)[:, 1]
                elseif field == "GT"
                    # Extract from GT (genotypes) field
                    allele_frequencies = sum(genotype_calls .== j, dims = 2)[:, 1] / ploidy
                else
                    throw(
                        ArgumentError(
                            "Unrecognized genotyped field: `" * field * "`. Please select `AF`, `AD` or `GT`.",
                        ),
                    )
                end
                # Set missing allele frequencies
                if length(idx_missing) > 0
                    allele_frequencies[idx_missing] .= missing
                end
                # Insert allele frequencies
                genomes.allele_frequencies[:, i] = allele_frequencies
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    close(file)
    # Checks
    unique_loci_alleles::Vector{String} = unique(genomes.loci_alleles)
    duplicated_loci_alleles::Vector{String} = []
    for locus_allele in unique_loci_alleles
        if sum(genomes.loci_alleles .== locus_allele) > 1
            push!(duplicated_loci_alleles, locus_allele)
        end
    end
    if length(duplicated_loci_alleles) > 0
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
