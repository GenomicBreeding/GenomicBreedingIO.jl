function vcf_count_loci(fname::String, verbose::Bool)::Tuple{Int64,Int64}
    if verbose
        println("Counting loci...")
    end
    n_loci::Int64 = 0
    total_lines::Int64 = 0
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
    for raw_line in eachline(file)
        # raw_line = readline(file)
        total_lines += 1
        if (raw_line[1] != '#') && (total_lines > 1) && (length(raw_line) > 0)
            n_loci += 1
        end
    end
    close(file)
    if verbose
        println(string("Total number of lines: ", total_lines))
        println(string("Total number of loci: ", n_loci))
    end
    total_lines, n_loci
end

function vcf_chunkify(
    fname::String,
    verbose::Bool,
    n_loci::Int64,
    line_counter::Int64,
)::Tuple{Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}}
    if verbose
        println("Counting threads and identifying the file positions and loci indexes per thread...")
    end
    n_threads = Threads.nthreads()
    n_lines_per_thread = Int(round(n_loci / n_threads))
    idx_loci_per_thread_ini = vcat([0], cumsum(repeat([n_lines_per_thread], n_threads - 1)))
    idx_loci_per_thread_fin = idx_loci_per_thread_ini .+ (n_lines_per_thread - 1)
    idx_loci_per_thread_fin[end] = n_loci
    file_pos_per_thread_ini = []
    file_pos_per_thread_fin = []
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
    counter::Int64 = 0
    previous_pos = position(file)
    for raw_line in eachline(file)
        # raw_line = readline(file)
        if (raw_line[1] != '#') && (line_counter > 1) && (length(raw_line) > 0)
            counter += 1
            if counter == 1
                append!(file_pos_per_thread_ini, previous_pos)
            end
            if counter == n_lines_per_thread
                append!(file_pos_per_thread_fin, position(file))
                counter = 0
            end
        end
        previous_pos = position(file)
    end
    if length(file_pos_per_thread_fin) == (n_threads - 1)
        try
            seekend(file)
        catch
            while !eof(file)
                readline(file)
            end
        end
        append!(file_pos_per_thread_fin, position(file))
    end
    # file_pos_per_thread_fin[end] = position(file)
    close(file)
    if verbose
        println(string("Total number of threads: ", n_threads))
    end
    (idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin)
end

function vcf_extract_entries_and_formats(fname::String, verbose::Bool)::Tuple{Vector{String},Vector{String}}
    # Column index of the start of genotype data
    IDX::Int64 = 10
    if verbose
        println("Extracting the names of the entries...")
    end
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
    header_line = ""
    for raw_line in eachline(file)
        if !isnothing(match(r"^#CHR", raw_line))
            header_line = raw_line
            break
        end
    end
    close(file)
    header::Vector{String} = split(header_line, "\t")
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
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
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
    if verbose
        println(string("Total number of entries: ", length(entries)))
    end
    (entries, format_lines)
end

function vcf_extract_info(
    fname::String,
    verbose::Bool,
    field::String,
    format_lines::Vector{String},
)::Tuple{String,Int64,Int64}
    # Column index of the start of genotype data
    IDX::Int64 = 10
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
    if verbose
        println("Identifying data field...")
    end
    format_details = split(format_lines[idx[1]], ",")
    field = split(format_details[.!isnothing.(match.(r"ID=", format_details))][1], "=")[end]
    n_alleles = 0
    ploidy = 0 # only for "GT" field - not used for "AD" and "AF" fields
    if field == "GT"
        # Computationally expensive allele counting for field=="GT":
        file =
            if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
                GZip.open(fname, "r")
            else
                open(fname, "r")
            end
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
        if verbose
            println("Field identified: \"GT\"")
            println(string("Ploidy: ", ploidy))
        end
    else
        # Easy n_alleles extraction for "AF" and "AD" fields
        n_alleles = try
            number = split(format_details[.!isnothing.(match.(r"Number=", format_details))][1], "=")[end]
            if number == "R"
                2
            else
                parse(Int64, number)
            end
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
        if verbose
            println("Field identified: " * field)
            println(string("Number of alleles per locus: ", n_alleles))
        end
    end
    (field, n_alleles, ploidy)
end

function vcf_instantiate_output(
    fname::String,
    verbose::Bool,
    entries::Vector{String},
    n_loci::Int64,
    n_alleles::Int64,
)::Genomes
    # Define the expected dimensions of the Genomes struct
    if verbose
        println("Initialising the output Genomes struct...")
    end
    n = length(entries)
    p = n_loci * (n_alleles - 1)
    # Instatiate the output struct
    genomes = Genomes(n = n, p = p)
    genomes.entries = entries
    genomes.populations .= "unknown"
    genomes.mask .= true
    if verbose
        println(string("Number of entries: ", n))
        println(string("Number of loci-alleles: ", p))
    end
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
    genomes
end

function vcf_parse_coordinates(
    line::Vector{String},
    line_counter::Int64,
    field::String,
    entries::Vector{String},
)::Union{Nothing,Tuple{Int64,String,Int64,String,Vector{String},Vector{String}}}
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
        return nothing
        if verbose
            println("Field `" * field * "` absent in line " * string(line_counter) * " of vcf file: `" * fname * "`.")
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
    refalts = if alt != ["."]
        vcat([ref], alt)
    else
        [ref]
    end
    (idx_field[1], chrom, pos, ref, alt, refalts)
end

function vcf_extract_allele_freqs!(
    genomes::Genomes,
    pb::Union{Nothing,Progress},
    i::Int64,
    fname::String,
    verbose::Bool,
    line_counter::Int64,
    field::String,
    idx_field::Int64,
    line::Vector{String},
    chrom::String,
    pos::Int64,
    ref::String,
    alt::Vector{String},
    refalts::Vector{String},
)::Tuple{Genomes,Union{Nothing,Progress},Int64}
    # Column index of the start of genotype data
    IDX::Int64 = 10
    if field == "AF"
        afreqs = try
            parse.(Float64, stack([split(split(x, ":")[idx_field], ",") for x in line[IDX:end]], dims = 1))
        catch
            throw(
                ErrorException(
                    "Cannot parse the `" *
                    field *
                    "` field (index=" *
                    string(idx_field) *
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
            parse.(Float64, stack([split(split(x, ":")[idx_field], ",") for x in line[IDX:end]], dims = 1))
        catch
            throw(
                ErrorException(
                    "Cannot parse the `" *
                    field *
                    "` field (index=" *
                    string(idx_field) *
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
            genotype_calls_tmp = vcat(split.(split(split(x, ":")[idx_field], "/"), "|")...)
            # Convert missing GT (i.e. '.') into zeroes
            replace!(genotype_calls_tmp, "." => "0")
            genotype_calls[k, :] = parse.(Int64, genotype_calls_tmp)
        end
        idx_missing = findall(sum(genotype_calls, dims = 2)[:, 1] .== 0.0)
    end
    for j in eachindex(refalts)
        # j = 2
        # Skip the reference allele
        if (j == 1) && (length(refalts) > 1)
            continue
        end
        allele = refalts[j]
        i += 1
        genomes.loci_alleles[i] = if length(refalts) > 1
            join([chrom, pos, join(refalts, "|"), allele], "\t")
        else
            join([chrom, pos, join(repeat(refalts, 2), "|"), allele], "\t")
        end
        allele_frequencies::Vector{Union{Missing,Float64}} = if field == "AF"
            # Extract from AF (allele frequencies) field
            afreqs[:, j]
        elseif field == "AD"
            # Extract from AD (allele depths) field
            depths[:, j] ./ sum(depths, dims = 2)[:, 1]
        elseif field == "GT"
            # Extract from GT (genotypes) field
            sum(genotype_calls .== j, dims = 2)[:, 1] / ploidy
        else
            throw(ArgumentError("Unrecognized genotyped field: `" * field * "`. Please select `AF`, `AD` or `GT`."))
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
    (genomes, pb, i)
end

"""
    readvcf(;fname::String, field::String = "any", verbose::Bool = false)::Genomes

Load Genomes struct from vcf file

# Examples
```jldoctest; setup = :(using GBCore, GBIO)
julia> genomes = GBCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> fname_gz = writevcf(genomes, gzip=true);

julia> genomes_reloaded = readvcf(fname=fname);

julia> genomes_reloaded_gz = readvcf(fname=fname_gz);

julia> genomes.entries == genomes_reloaded.entries == genomes_reloaded_gz.entries
true

julia> dimensions(genomes) == dimensions(genomes_reloaded) == dimensions(genomes_reloaded_gz)
true

julia> ismissing.(genomes.allele_frequencies) == ismissing.(genomes_reloaded.allele_frequencies) == ismissing.(genomes_reloaded_gz.allele_frequencies)
true
```
"""
function readvcf(; fname::String, field::String = "any", verbose::Bool = false)::Genomes
    # genomes = GBCore.simulategenomes(n=10, sparsity=0.1); fname = writevcf(genomes, gzip=true); field = "any"; verbose = true;
    # genomes = simulategenomes(n_alleles=3, sparsity=0.1); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 4) ./ 4; fname = writevcf(genomes, ploidy=4); field = "GT"; verbose = true;
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    gzip = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        true
    else
        false
    end
    # Count the number of lines in the file which are not header lines or comments
    total_lines, n_loci = vcf_count_loci(fname, verbose)
    # Find location of each file chunk for multi-threaded parsing
    (idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin) =
        vcf_chunkify(fname, verbose, n_loci, total_lines)
    # Extract the names of the entries or samples
    entries, format_lines = vcf_extract_entries_and_formats(fname, verbose)
    # Extract info
    field, n_alleles, ploidy = vcf_extract_info(fname, verbose, field, format_lines)
    # Instantiate the output Genomes struct
    genomes = vcf_instantiate_output(fname, verbose, entries, n_loci, n_alleles)
    n, p = size(genomes.allele_frequencies)
    # Read the file line by line
    pb = if verbose
        ProgressMeter.Progress(p; desc = "Loading genotypes from vcf file: ")
    else
        nothing
    end
    files = []
    for i in eachindex(file_pos_per_thread_ini)
        file =
            if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
                GZip.open(fname, "r")
            else
                open(fname, "r")
            end
        push!(files, file)
    end
    thread_lock::ReentrantLock = ReentrantLock()
    Threads.@threads for idx in eachindex(file_pos_per_thread_ini)
        # for idx in eachindex(file_pos_per_thread_ini)
        # idx = 2
        # println(string("idx = ", idx))
        allele_frequencies::Vector{Union{Missing,Float64}} = fill(missing, n)
        ini = file_pos_per_thread_ini[idx]
        fin = file_pos_per_thread_fin[idx]
        line::Vector{String} = repeat([""], n + 9)
        file = files[idx]
        seek(file, ini)
        i = idx_loci_per_thread_ini[idx]
        line_counter = (i + (total_lines - n_loci)) - 1
        while (position(file) <= fin) && !eof(file)
            raw_line = readline(file)
            line_counter += 1
            # If we somehow end the end of the file or encounter an empty line
            if (length(raw_line) == 0) && (position(file) == fin)
                break # end of file
            elseif (length(raw_line) == 0) && (position(file) < fin)
                continue # empty line in the body of the file
            end
            # Skip commented out lines including the first 2 header
            if raw_line[1] == '#'
                continue
            end
            line .= split(raw_line, "\t")
            info = vcf_parse_coordinates(line, line_counter, field, entries)
            if isnothing(info)
                continue
            end
            idx_field = info[1]
            chrom = info[2]
            pos = info[3]
            ref = info[4]
            alt = info[5]
            refalts = info[6]
            @lock thread_lock genomes, pb, i = vcf_extract_allele_freqs!(
                genomes,
                pb,
                i,
                fname,
                verbose,
                line_counter,
                field,
                idx_field,
                line,
                chrom,
                pos,
                ref,
                alt,
                refalts,
            )
        end
        close(file)
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    # Remove missing loci-alleles
    idx = findall(genomes.loci_alleles .!= "")
    if verbose && (length(idx) < length(genomes.loci_alleles))
        println(string("Removing ", length(genomes.loci_alleles) - length(idx), " empty loci-alleles..."))
    end
    genomes.loci_alleles = genomes.loci_alleles[idx]
    genomes.allele_frequencies = genomes.allele_frequencies[:, idx]
    genomes.mask = genomes.mask[:, idx]
    # Checks
    if verbose
        println("Output checks...")
    end
    unique_loci_alleles::Vector{String} = unique(genomes.loci_alleles)
    duplicated_loci_alleles::Vector{String} = []
    if length(unique_loci_alleles) < length(genomes.loci_alleles)
        for locus_allele in unique_loci_alleles
            if sum(genomes.loci_alleles .== locus_allele) > 1
                push!(duplicated_loci_alleles, locus_allele)
            end
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
