"""
    vcfcountlocialleles(fname::String; verbose::Bool = false)::Tuple{Int64,Int64}

Count the number of loci and total lines in a VCF file.

# Arguments
- `fname::String`: Path to the VCF file. Can be either a plain text VCF file or a gzipped VCF file (with extensions .vcf.gz or .vcf.bgz)
- `verbose::Bool`: If true, prints progress messages and results to stdout. Defaults to false.

# Returns
- `Tuple{Int64,Int64}`: A tuple containing:
    - First element: Total number of lines in the file (including headers)
    - Second element: Number of data lines (variants/loci) excluding header lines

# Description
Reads through a VCF (Variant Call Format) file and counts:
1. Total lines in the file (including headers)
2. Number of data lines (variants/loci) that don't start with '#'

The function automatically detects and handles different file formats:
- Plain text VCF files (.vcf)
- Gzipped VCF files (.vcf.gz)
- BGZipped VCF files (.vcf.bgz)

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> fname_gz = writevcf(genomes, gzip=true);

julia> n_1, p_1, l_1 = vcfcountlocialleles(fname);

julia> n_2, p_2, l_2 = vcfcountlocialleles(fname_gz);

julia> n_1 == n_2 == 10_009
true

julia> p_1 == p_2 == 10_000
true

julia> l_1 == l_2 == 10_000
true

julia> rm.([fname, fname_gz]);
```
"""
function vcfcountlocialleles(fname::String; verbose::Bool = false)::Tuple{Int64,Int64,Int64}
    if verbose
        println("Counting loci...")
    end
    n_loci::Int64 = 0
    n_alt_alleles::Int64 = 0
    total_lines::Int64 = 0
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
    for raw_line in eachline(file)
        # raw_line = readline(file)
        total_lines += 1
        if (raw_line[1] != '#') && (length(raw_line) > 0)
            n_loci += 1
            n_alt_alleles += length(split(split(raw_line, "\t")[5], ",")) # alt alleles
        end
    end
    close(file)
    if verbose
        println(string("Total number of lines: ", total_lines))
        println(string("Total number of loci: ", n_loci))
    end
    total_lines, n_loci, n_alt_alleles
end

"""
    vcfchunkify(fname::String; n_loci::Int64, verbose::Bool = false)::Tuple{Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}}

Divide a VCF file into chunks for parallel processing.

# Arguments
- `fname::String`: Path to the VCF file (can be .vcf, .vcf.gz, or .vcf.bgz)
- `n_loci::Int64`: Total number of loci in the VCF file
- `verbose::Bool=false`: If true, prints progress information

# Returns
A tuple containing four Vector{Int64} arrays:
1. Starting loci indices for each thread
2. Ending loci indices for each thread
3. Starting file positions for each thread
4. Ending file positions for each thread

# Details
- Automatically detects if the input file is gzipped
- Divides the workload evenly across available threads
- Skips header lines (starting with '#')
- Handles both regular and gzipped VCF files

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> _, n_loci, n_alt_alleles = vcfcountlocialleles(fname);

julia> idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin = vcfchunkify(fname, n_loci=n_loci);

julia> length(idx_loci_per_thread_ini) == length(idx_loci_per_thread_fin) == length(file_pos_per_thread_ini) == length(file_pos_per_thread_fin)
true

julia> (idx_loci_per_thread_ini[1] == 0) && (sum(idx_loci_per_thread_ini .== 0) == 1)
true

julia> (idx_loci_per_thread_fin[end] == n_loci) && (sum(idx_loci_per_thread_fin .== 0) == 0)
true

julia> (sum(file_pos_per_thread_ini .== 0) == 0) && (sum(file_pos_per_thread_fin .== 0) == 0)
true

julia> rm(fname);
```
"""
function vcfchunkify(
    fname::String;
    n_loci::Int64,
    verbose::Bool = false,
)::Tuple{Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}}
    if verbose
        println("Counting threads and identifying the file positions and loci indexes per thread...")
    end
    n_threads = Threads.nthreads()
    n_lines_per_thread = Int(ceil(n_loci / n_threads))
    idx_loci_per_thread_ini = vcat([0], cumsum(repeat([n_lines_per_thread], n_threads - 1)))
    idx_loci_per_thread_fin = idx_loci_per_thread_ini .+ (n_lines_per_thread - 1)
    idx_loci_per_thread_fin[end] = n_loci
    file_pos_per_thread_ini::Vector{Int64} = zeros(n_threads)
    file_pos_per_thread_fin::Vector{Int64} = zeros(n_threads)
    file = if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
        GZip.open(fname, "r")
    else
        open(fname, "r")
    end
    i = 1
    counter::Int64 = 0
    previous_pos = position(file)
    for raw_line in eachline(file)
        # raw_line = readline(file)
        if (raw_line[1] != '#') && (length(raw_line) > 0)
            counter += 1
            if counter == 1
                file_pos_per_thread_ini[i] = previous_pos
            end
            if counter == n_lines_per_thread
                file_pos_per_thread_fin[i] = position(file)
                i += 1
                counter = 0
            end
        end
        previous_pos = position(file)
    end
    if file_pos_per_thread_fin[end] == 0
        try
            seekend(file)
        catch
            while !eof(file)
                readline(file)
            end
        end
        file_pos_per_thread_fin[i] = position(file)
    end
    close(file)
    if verbose
        println(string("Total number of threads: ", n_threads))
    end
    (idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin)
end

"""
    vcfextractentriesandformats(fname::String; verbose::Bool = false)::Tuple{Vector{String},Vector{String}}

Extract sample entries and format definitions from a VCF file.

# Arguments
- `fname::String`: Path to the VCF file (can be gzipped with extensions .vcf.gz or .vcf.bgz)
- `verbose::Bool=false`: If true, prints progress information to stdout

# Returns
A tuple containing:
1. `Vector{String}`: List of sample names from the VCF header
2. `Vector{String}`: List of FORMAT field definitions from the VCF metadata

# Description
Reads a VCF file and extracts two key pieces of information:
1. Sample names from the header line (columns after FORMAT field)
2. FORMAT field definitions from metadata lines starting with "##FORMAT"

The function validates the presence and correct order of standard VCF columns:
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT

# Throws
- `ArgumentError`: If VCF has fewer than expected columns or column names don't match VCF format

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> entries, format_lines = vcfextractentriesandformats(fname);

julia> entries
10-element Vector{String}:
 "entry_01"
 "entry_02"
 "entry_03"
 "entry_04"
 "entry_05"
 "entry_06"
 "entry_07"
 "entry_08"
 "entry_09"
 "entry_10"

julia> format_lines
3-element Vector{String}:
 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
 "##FORMAT=<ID=AD,Number=2,Type=Float,Description=\\"Allele Depth\\">"
 "##FORMAT=<ID=AF,Number=2,Type=Float,Description=\\"Allele Frequency\\">"

julia> rm(fname);
```
"""
function vcfextractentriesandformats(fname::String; verbose::Bool = false)::Tuple{Vector{String},Vector{String}}
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
    # Output
    (entries, format_lines)
end


"""
    vcfextractinfo(fname::String; format_lines::Vector{String}, field::String="any", verbose::Bool=false)::Tuple{String,Int64,Int64}

Extract information about genotype fields from a VCF file.

# Arguments
- `fname::String`: Path to the VCF file (can be gzipped)
- `format_lines::Vector{String}`: Vector containing FORMAT lines from the VCF header
- `field::String="any"`: Specific field to extract ("GT", "AD", "AF", or "any")
- `verbose::Bool=false`: If true, prints progress information

# Returns
A tuple containing:
1. `field::String`: The identified genotype field
2. `n_alleles::Int64`: Maximum number of alleles per locus
3. `ploidy::Int64`: Ploidy level (only meaningful for GT field; set to typemax(Int64) for AD and AF fields)

# Details
- If `field` is "any", searches for fields in priority order: AF > AD > GT
- For GT field, scans entire file to determine maximum number of alleles and ploidy
- For AF and AD fields, extracts allele count from format header
- Supports both gzipped (.gz, .bgz) and uncompressed VCF files

# Throws
- `ArgumentError`: If specified field is not found in the VCF file
- `ErrorException`: If unable to parse number of alleles from format header

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> _, format_lines = vcfextractentriesandformats(fname);

julia> field, n_alleles, ploidy = vcfextractinfo(fname, format_lines=format_lines);

julia> (field == "AF") && (n_alleles == 2) && (ploidy == typemax(Int64))
true

julia> rm(fname);
```
"""
function vcfextractinfo(
    fname::String;
    format_lines::Vector{String},
    field::String = "any",
    verbose::Bool = false,
)::Tuple{String,Int64,Int64}
    # Column index of the start of genotype data
    IDX::Int64 = 10
    # Choose field where priority order starts with AF with the highest priority followed by AD, and finally GT
    idx = if field == "any"
        idx = findall(.!isnothing.(match.(r"ID=AF", format_lines)))
        if length(idx) == 0
            idx = findall(.!isnothing.(match.(r"ID=AD", format_lines)))
        end
        if length(idx) == 0
            idx = findall(.!isnothing.(match.(r"ID=GT", format_lines)))
        end
        if length(idx) == 0
            throw(ArgumentError("The input vcf file: `" * fname * "` does not have `AF`, `AD`, or `GT` genotype fields."))
        end
        idx
    else
        idx = findall(.!isnothing.(match.(Regex(field), format_lines)))
        if length(idx) == 0
            throw(ArgumentError("The input vcf file: `" * fname * "` does not have `" * field * "` field."))
        end
        idx
    end
    # Define the field and the maximum number of alleles per locus
    if verbose
        println("Identifying data field...")
    end
    format_details = split(format_lines[idx[1]], ",")
    field = split(format_details[.!isnothing.(match.(r"ID=", format_details))][1], "=")[end]
    n_alleles, ploidy = if field == "GT"
        # Computationally expensive allele counting for field=="GT":
        file =
            if (split(fname, ".")[(end-1):end] == ["vcf", "gz"]) || (split(fname, ".")[(end-1):end] == ["vcf", "bgz"])
                GZip.open(fname, "r")
            else
                open(fname, "r")
            end
        n_alleles = 0
        ploidy = 0
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
            # Ploidy only for "GT" field - not used for "AD" and "AF" fields
            if ploidy < ψ
                ploidy = ψ
            end
        end
        close(file)
        if verbose
            println("Field identified: \"GT\"")
            println(string("Ploidy: ", ploidy))
        end
        n_alleles, ploidy
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
            println(
                string(
                    "Number of alleles per locus according to FORMAT lines (fret not though as number of alleles is counted per locus and may not correspond to this): ",
                    n_alleles,
                ),
            )
        end
        n_alleles, typemax(Int64)
    end
    # Output
    (field, n_alleles, ploidy)
end

"""
    vcfinstantiateoutput(fname::String; entries::Vector{String}, n_alt_alleles::Int64, verbose::Bool = false)::Genomes

Create and initialize a Genomes struct from VCF file parameters.

# Arguments
- `fname::String`: Name of the VCF file being processed
- `entries::Vector{String}`: Vector containing entry identifiers
- `n_alt_alleles::Int64`: Total number of alternative alleles across all loci
- `verbose::Bool=false`: If true, prints progress information

# Returns
- `Genomes`: An initialized Genomes struct with:
  - dimensions n × p where n is number of entries and p = n_alt_alleles
  - entry names assigned from input entries
  - populations set to "unknown"
  - mask set to true
  
# Throws
- `ErrorException`: If duplicate entries are found in the VCF file

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> entries, format_lines = vcfextractentriesandformats(fname);

julia> _, _, n_alt_alleles = vcfcountlocialleles(fname);

julia> genomes_instantiated = vcfinstantiateoutput(fname, entries=entries, n_alt_alleles=n_alt_alleles);

julia> size(genomes_instantiated.allele_frequencies)
(10, 10000)

julia> genomes_instantiated.entries == entries
true

julia> rm(fname);
```
"""
function vcfinstantiateoutput(
    fname::String;
    entries::Vector{String},
    n_alt_alleles::Int64,
    verbose::Bool = false,
)::Genomes
    # Define the expected dimensions of the Genomes struct
    if verbose
        println("Initialising the output Genomes struct...")
    end
    n = length(entries)
    p = n_alt_alleles # equivalent to n_loci * (n_alleles - 1) if n_alleles are the same across all loci
    # Instantiate the output struct
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

"""
    vcfparsecoordinates(; line::Vector{String}, line_counter::Int64, field::String)::Union{Nothing,Tuple{Int64,String,Int64,Vector{String}}}

Parse coordinates and allele information from a VCF file line.

# Arguments
- `line::Vector{String}`: A vector containing the split line from VCF file
- `line_counter::Int64`: Current line number being processed in the VCF file
- `field::String`: The field name to extract allele frequencies from

# Returns
- `Nothing`: If the specified field is not found in the line
- `Tuple{Int64,String,Int64,Vector{String}}`: A tuple containing:
    - Field index
    - Chromosome name
    - Position
    - Combined reference and alternative alleles

# Throws
- `ErrorException`: If the position field cannot be parsed as an integer

# Note
The function validates the line format and extracts genomic coordinates and allele information
from a VCF file line. It handles missing alternative alleles (denoted by ".") and performs
necessary type conversions.

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> entries, format_lines = vcfextractentriesandformats(fname);

julia> field, _, _ = vcfextractinfo(fname, format_lines=format_lines);

julia> file = open(fname, "r"); line::Vector{String} = split([readline(file) for i in 1:10][end], "\t"); close(file);

julia> idx_field, chrom, pos, refalts = vcfparsecoordinates(line=line, line_counter=10, field=field);

julia> (idx_field == 3) && (chrom == line[1]) && (pos == parse(Int64, line[2])) && (refalts == line[4:5])
true

julia> rm(fname);
```
"""
function vcfparsecoordinates(;
    line::Vector{String},
    line_counter::Int64,
    field::String,
)::Union{Nothing,Tuple{Int64,String,Int64,Vector{String}}}
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
    (idx_field[1], chrom, pos, refalts)
end

"""
    vcfextractallelefreqs!(genomes::Genomes, pb::Union{Nothing,Progress}, i::Vector{Int64}; 
                          fname::String, line::Vector{String}, line_counter::Int64, 
                          field::String, min_depth::Int64=10, max_depth::Int64=100, 
                          verbose::Bool=false)

Extract allele frequencies from VCF file data and update a Genomes object.

# Arguments
- `genomes::Genomes`: Object to store genomic data
- `pb::Union{Nothing,Progress}`: Progress bar object or nothing
- `i::Vector{Int64}`: Single-element vector containing current locus-allele index
- `fname::String`: Name of VCF file being processed
- `line::Vector{String}`: Current line from VCF file split into fields
- `line_counter::Int64`: Current line number in VCF file
- `field::String`: Type of field to extract ("AF", "AD", or "GT")
- `min_depth::Int64=10`: Minimum read depth threshold for AD field
- `max_depth::Int64=100`: Maximum read depth threshold for AD field
- `verbose::Bool=false`: Whether to display progress updates

# Returns
Nothing; Updates the input parameters in place:
- `genomes`: Updated with new allele frequencies and loci information
- `pb`: Advanced if verbose=true
- `i`: Index incremented based on processed alleles

# Description
Processes VCF data to extract allele frequencies of the alternative allele/s using one of three methods:
- AF field: Direct frequency values from VCF
- AD field: Calculated from read depths (filtered by min_depth and max_depth)
- GT field: Calculated from genotype calls

Updates Genomes object with:
- Loci-allele identifiers (chromosome, position, alleles)
- Allele frequencies for each sample

# Throws
- `ArgumentError`: If field parameter is not "AF", "AD", or "GT"
- `ErrorException`: If unable to parse AF or AD fields from VCF

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

julia> fname = writevcf(genomes);

julia> _, _, n_alt_alleles = vcfcountlocialleles(fname);

julia> entries, format_lines = vcfextractentriesandformats(fname);

julia> field, n_alleles, _ = vcfextractinfo(fname, format_lines=format_lines);

julia> genomes_instantiated = vcfinstantiateoutput(fname, entries=entries, n_alt_alleles=n_alt_alleles);

julia> sum(ismissing.(genomes_instantiated.allele_frequencies[:, 1])) == length(entries)
true

julia> file = open(fname, "r"); line::Vector{String} = split([readline(file) for i in 1:10][end], "\t"); close(file);

julia> vcfextractallelefreqs!(genomes_instantiated, nothing, [0], fname=fname, line=line, line_counter=10, field=field);

julia> sum(ismissing.(genomes_instantiated.allele_frequencies[:, 1])) == 0
true

julia> rm(fname);
```
"""
function vcfextractallelefreqs!(
    genomes::Genomes,
    pb::Union{Nothing,Progress},
    i::Vector{Int64};
    fname::String,
    line::Vector{String},
    line_counter::Int64,
    field::String,
    min_depth::Int64 = 10,
    max_depth::Int64 = 100,
    verbose::Bool = false,
)
    # Parse coordinates
    info = vcfparsecoordinates(line = line, line_counter = line_counter, field = field)
    if isnothing(info)
        # Skip if AF, AD, and GT fields are all empty
        return nothing
    end
    idx_field = info[1]
    chrom = info[2]
    pos = info[3]
    refalts = info[4]
    # Column index of the start of genotype data
    IDX::Int64 = 10
    (afreqs, depths, genotype_calls, idx_missing) = if field == "AF"
        afreqs = try
            parse.(
                Float64,
                stack(
                    [[af == "." ? "0.0" : af for af in split(split(x, ":")[idx_field], ",")] for x in line[IDX:end]],
                    dims = 1,
                ),
            )
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
        (afreqs, [], [], idx_missing)
    elseif field == "AD"
        # Parse allele depths (AD field) from the VCF file
        depths = try
            # Extract AD field values for each sample, replacing missing values (".") with "0"
            ads = [[ad == "." ? "0" : string(ad) for ad in split(split(x, ":")[idx_field], ",")] for x in line[IDX:end]]
            # Determine the maximum number of alleles across all samples
            a_max = maximum([length(x) for x in ads])
            # Pad allele depth arrays with zeros to ensure consistent dimensions across samples
            ads = [length(x) < a_max ? vcat(x, repeat(["0"], a_max - length(x))) : x for x in ads]
            # Convert allele depth strings to Float64 and stack them into a matrix
            parse.(Float64, stack(ads, dims = 1))
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
        # Set depth beyond the min and max depth to zero
        depths[(depths.<min_depth).||(depths.>max_depth)] .= 0.0
        idx_missing = findall(sum(depths, dims = 2)[:, 1] .== 0.0)
        ([], depths, [], idx_missing)
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
        ([], [], genotype_calls, idx_missing)
    end
    for j in eachindex(refalts)
        # j = 2
        # Skip the reference allele (i.e. the first allele)
        if (j == 1) && (length(refalts) > 1)
            continue
        end
        allele = refalts[j]
        i[1] += 1
        genomes.loci_alleles[i[1]] = if length(refalts) > 1
            # Biallelic or multi-allelic instances (usual/expected case)
            join([chrom, pos, join(refalts, "|"), allele], "\t")
        else
            # Mono-allelic instances
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
        genomes.allele_frequencies[:, i[1]] = allele_frequencies
        if verbose
            ProgressMeter.next!(pb)
        end
    end
end

"""
    readvcf(; fname::String, field::String = "any", min_depth::Int64 = 5, max_depth::Int64 = 100, verbose::Bool = false)::Genomes

Read genetic data from a VCF (Variant Call Format) file into a Genomes struct.

# Arguments
- `fname::String`: Path to the VCF file. Can be gzipped (.vcf.gz or .vcf.bgz) or uncompressed (.vcf)
- `field::String="any"`: Which FORMAT field to extract from VCF. Default "any" tries to automatically detect genotype field
- `min_depth::Int64=5`: Minimum read depth threshold for AD (Allele Depth) field
- `max_depth::Int64=100`: Maximum read depth threshold for AD field
- `verbose::Bool=false`: Whether to print progress and debug information

# Returns
- `Genomes`: A Genomes struct containing the loaded genetic data with fields:
  - `allele_frequencies`: Matrix of allele frequencies
  - `loci_alleles`: Vector of locus-allele combination strings
  - `mask`: Boolean matrix indicating missing data
  - `samples`: Vector of sample names

# Details
Reads VCF files in parallel using multiple threads. Handles multi-allelic variants and different ploidies. 
Field priority (when field="any"):
1. AF (Allele Frequency)  
2. AD (Allele Depth)  
3. GT (Genotype)

Performs various checks on the input data including:
- File existence
- No duplicate loci-allele combinations 
- Consistent dimensions in output struct

# Throws
- `ErrorException`: If file doesn't exist, has duplicates, or output dimensions are invalid

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);

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
function readvcf(;
    fname::String,
    field::String = "any",
    min_depth::Int64 = 10,
    max_depth::Int64 = 100,
    verbose::Bool = false,
)::Genomes
    # genomes = simulategenomes(n_alleles=3, sparsity=0.1); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 4) ./ 4; fname = writevcf(genomes, ploidy=4); field = "GT"; verbose = true;
    # Check input arguments
    if !isfile(fname)
        throw(ErrorException("The file: " * fname * " does not exist."))
    end
    # Count the number of lines in the file which are not header lines or comments
    total_lines, n_loci, n_alt_alleles = vcfcountlocialleles(fname, verbose = verbose)
    # Find location of each file chunk for multi-threaded parsing
    (idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin) =
        vcfchunkify(fname, n_loci = n_loci, verbose = verbose)
    # Extract the names of the entries or samples
    entries, format_lines = vcfextractentriesandformats(fname, verbose = verbose)
    # Extract info
    field, n_alleles, ploidy = vcfextractinfo(fname, field = field, format_lines = format_lines, verbose = verbose)
    # Instantiate the output Genomes struct
    genomes = vcfinstantiateoutput(fname, entries = entries, n_alt_alleles = n_alt_alleles, verbose = verbose)
    n, p = size(genomes.allele_frequencies)
    # Instantiate the progress meter
    pb = if verbose
        ProgressMeter.Progress(p; desc = "Loading genotypes from vcf file: ")
    else
        nothing
    end
    # Create a vector of open file streams, one for each thread
    files = []
    for _ in eachindex(file_pos_per_thread_ini)
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
        # idx = 1
        # println(string("idx = ", idx))
        ini = file_pos_per_thread_ini[idx]
        fin = file_pos_per_thread_fin[idx]
        file = files[idx]
        seek(file, ini)
        i::Vector{Int64} = [idx_loci_per_thread_ini[idx]]
        line_counter = (i[1] + (total_lines - n_loci)) - 1
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
            line::Vector{String} = split(raw_line, "\t")
            # println(line)
            if (n + 9) != length(line)
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
            @lock thread_lock vcfextractallelefreqs!(
                genomes,
                pb,
                i,
                fname = fname,
                line = line,
                line_counter = line_counter,
                field = field,
                min_depth = min_depth,
                max_depth = max_depth,
                verbose = verbose,
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
    writevcf(genomes::Genomes; fname::Union{Missing,String} = missing, ploidy::Int64 = 0, 
             max_depth::Int64 = 100, n_decimal_places::Int64 = 4, gzip::Bool = false)::String

Write genomic data to a Variant Call Format (VCF) file.

# Arguments
- `genomes::Genomes`: A Genomes object containing the genetic data to be written.
- `fname::Union{Missing,String} = missing`: Output filename. If missing, generates a default name with timestamp.
- `ploidy::Int64 = 0`: The ploidy level of the organisms (e.g., 2 for diploid).
- `max_depth::Int64 = 100`: Maximum depth for allele depth calculation.
- `n_decimal_places::Int64 = 4`: Number of decimal places for rounding allele frequencies.
- `gzip::Bool = false`: Whether to compress the output file using gzip.

# Returns
- `String`: The name of the created VCF file.

# Description
Creates a VCF v4.2 format file containing genomic variants data. The function processes
allele frequencies and depths, calculates genotypes based on ploidy, and formats the
data according to VCF specifications. The output includes:
- Standard VCF header information
- Sample information with FORMAT fields:
  * GT (Genotype)
  * AD (Allele Depth)
  * AF (Allele Frequency)

# Throws
- `DimensionMismatch`: If the input Genomes object has inconsistent dimensions
- `ErrorException`: If the output file already exists
- `ArgumentError`: If the file extension is not '.vcf' or if the specified directory doesn't exist

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO)
julia> genomes_1 = GenomicBreedingCore.simulategenomes(n=2, verbose=false);

julia> writevcf(genomes_1, fname="test_genomes_1.vcf")
"test_genomes_1.vcf"

julia> genomes_2 = GenomicBreedingCore.simulategenomes(n=2, n_alleles=3, verbose=false);

julia> genomes_2.allele_frequencies = round.(genomes_2.allele_frequencies .* 4) ./ 4;

julia> writevcf(genomes_2, fname="test_genomes_2.vcf", ploidy=4)
"test_genomes_2.vcf"

julia> genomes_3 = GenomicBreedingCore.simulategenomes(n=3, verbose=false);

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
    overwrite::Bool = false,
    verbose::Bool = false,
)::String
    # genomes = simulategenomes(n_alleles=3, sparsity=0.10); fname = missing; ploidy = 0; max_depth = 100; n_decimal_places = 4; gzip = true;
    # genomes = simulategenomes(n_alleles=3); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 2) ./ 2; fname = missing; ploidy = 2; max_depth = 100; n_decimal_places = 4; gzip = true;
    # genomes = simulategenomes(n_alleles=3); genomes.allele_frequencies = round.(genomes.allele_frequencies .* 4) ./ 4; fname = missing; ploidy = 4; max_depth = 100; n_decimal_places = 4; gzip = true;
    # Check input arguments
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted ☹."))
    end
    if ismissing(fname)
        fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".vcf")
    else
        if isfile(fname)
            if overwrite
                rm(fname)
            else
                throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
            end
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
        "##source=GenomicBreedingIO.jl-v1.0.0-DEV",
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
    if verbose
        pb = ProgressMeter.Progress(length(loci_ini_idx); desc = "Writing VCF file: ")
    end
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
        if verbose
            ProgressMeter.next!(pb)
        end
    end
    close(file)
    if verbose
        ProgressMeter.finish!(pb)
    end
    return fname
end
