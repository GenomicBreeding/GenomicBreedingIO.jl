var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GenomicBreedingIO","category":"page"},{"location":"#GenomicBreedingIO","page":"Home","title":"GenomicBreedingIO","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GenomicBreedingIO.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GenomicBreedingIO]","category":"page"},{"location":"#GenomicBreedingIO.isfuzzymatch-Tuple{String, String}","page":"Home","title":"GenomicBreedingIO.isfuzzymatch","text":"isfuzzymatch(a::String, b::String; threshold::Float64=0.3)::Bool\n\nDetermines if two strings approximately match each other using Levenshtein distance.\n\nThe function compares two strings and returns true if they are considered similar enough based on the Levenshtein edit distance and a threshold value. The threshold is applied as a fraction of the length of the shorter string.\n\nArguments\n\na::String: First string to compare\nb::String: Second string to compare\nthreshold::Float64=0.3: Maximum allowed edit distance as a fraction of the shorter string length\n\nReturns\n\nBool: true if the strings match within the threshold, false otherwise\n\nExamples\n\njulia> isfuzzymatch(\"populations\", \"populations\")\ntrue\n\njulia> isfuzzymatch(\"populations\", \"poplation\")\ntrue\n\njulia> isfuzzymatch(\"populations\", \"entry\")\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.levenshteindistance-Tuple{String, String}","page":"Home","title":"GenomicBreedingIO.levenshteindistance","text":"levenshteindistance(a::String, b::String)::Int64\n\nCalculate the Levenshtein distance (edit distance) between two strings.\n\nThe Levenshtein distance is a measure of the minimum number of single-character edits  (insertions, deletions, or substitutions) required to change one string into another.\n\nArguments\n\na::String: First input string\nb::String: Second input string\n\nReturns\n\nInt64: The minimum number of edits needed to transform string a into string b\n\nExamples\n\njulia> levenshteindistance(\"populations\", \"populations\")\n0\n\njulia> levenshteindistance(\"populations\", \"poplation\")\n2\n\njulia> levenshteindistance(\"populations\", \"entry\")\n3\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.readdelimited-Tuple{Type{GenomicBreedingCore.Genomes}}","page":"Home","title":"GenomicBreedingIO.readdelimited","text":"readdelimited(\n    type::Type{Genomes};\n    fname::String,\n    sep::String = \"\\t\",\n    parse_populations_from_entries::Union{Nothing,Function} = nothing,\n    verbose::Bool = false\n)::Genomes\n\nLoad genotype data from a delimited text file into a Genomes struct.\n\nArguments\n\ntype::Type{Genomes}: Type parameter (always Genomes)\nfname::String: Path to the input file\nsep::String: Delimiter character (default: tab)\nparse_populations_from_entries::Union{Nothing,Function}: Optional function to extract population names from entry names\nverbose::Bool: Whether to show progress bar during loading\n\nFile Format\n\nThe input file should be structured as follows:\n\nSupported extensions: .tsv, .csv, or .txt\nComments and headers start with '#'\nHeader format (2 lines where the second line is optional):\nColumn names: \"chrom,pos,allalleles,allele,entry1,entry_2,...\"\nPopulation names (optional): \"chrom,pos,allalleles,allele,pop1,pop_2,...\"\nData columns:\nchromosome identifier\nposition (integer)\nall alleles at locus (delimited by '|')\nspecific allele\n5+. allele frequencies for each entry (0.0-1.0 or missing/NA)\n\nReturns\n\nGenomes: A populated Genomes struct containing the loaded data\n\nThrows\n\nErrorException: If file doesn't exist or has invalid format\nArgumentError: If column names don't match expected format\nOverflowError: If allele frequencies are outside [0,1] range\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> genomes.entries = [string(genomes.populations[i], \"-\", genomes.entries[i]) for i in eachindex(genomes.populations)];\n\njulia> fname = writedelimited(genomes);\n\njulia> genomes_reloaded = readdelimited(Genomes, fname=fname);\n\njulia> genomes == genomes_reloaded\ntrue\n\njulia> fname = writedelimited(genomes, include_population_header=false);\n\njulia> genomes_reloaded = readdelimited(Genomes, fname=fname);\n\njulia> unique(genomes_reloaded.populations) == [\"Unknown_population\"]\ntrue\n\njulia> genomes_reloaded = readdelimited(Genomes, fname=fname, parse_populations_from_entries=x -> split(x, \"-\")[1]);\n\njulia> genomes == genomes_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.readdelimited-Tuple{Type{GenomicBreedingCore.Phenomes}}","page":"Home","title":"GenomicBreedingIO.readdelimited","text":"readdelimited(type::Type{Phenomes}; fname::String, sep::String = \"\\t\", verbose::Bool = false)::Phenomes\n\nLoad phenotypic data from a delimited text file into a Phenomes struct.\n\nArguments\n\ntype::Type{Phenomes}: Type parameter (must be Phenomes)\nfname::String: Path to the input file\nsep::String: Delimiter character (default: tab \"\\t\")\nverbose::Bool: Whether to show progress bar during loading (default: false)\n\nFile Format\n\nThe file should be a delimited text file with:\n\nHeader row containing column names\nFirst column: Entry identifiers\nSecond column: Population identifiers \nRemaining columns: Phenotypic trait values (numeric or missing)\n\nMissing values can be specified as \"missing\", \"NA\", \"na\", \"N/A\", \"n/a\" or empty string.\n\nReturns\n\nPhenomes: A Phenomes struct containing the loaded phenotypic data\n\nThrows\n\nErrorException: If file doesn't exist or has invalid format\nArgumentError: If required columns are missing or misnamed\nErrorException: If duplicate entries or traits are found\nErrorException: If numeric values cannot be parsed\n\nNotes\n\nComments starting with '#' are ignored\nEmpty lines are skipped\nMathematical operators (+,-,*,/,%) in trait names are replaced with underscores\nPerforms dimension checks on the loaded data\n\nExamples\n\njulia> phenomes = Phenomes(n=10, t=3); phenomes.entries = string.(\"entry_\", 1:10); phenomes.populations .= \"pop1\"; phenomes.traits = [\"A\", \"B\", \"C\"]; phenomes.phenotypes = rand(10,3); phenomes.phenotypes[1,1] = missing; phenomes.mask .= true;\n\njulia> fname = writedelimited(phenomes);\n\njulia> phenomes_reloaded = readdelimited(Phenomes, fname=fname);\n\njulia> phenomes == phenomes_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.readdelimited-Tuple{Type{GenomicBreedingCore.Trials}}","page":"Home","title":"GenomicBreedingIO.readdelimited","text":"readdelimited(type::Type{Trials}; fname::String, sep::String = \"\\t\", verbose::Bool = false)::Trials\n\nLoad a Trials struct from a string-delimited file.\n\nArguments\n\ntype::Type{Trials}: Type parameter (must be Trials)\nfname::String: Path to the input file\nsep::String = \"\\t\": Delimiter character (default is tab)\nverbose::Bool = false: Whether to display progress information\n\nRequired File Structure\n\nThe input file must contain the following 10 identifier columns:\n\nyears: Year identifiers\nseasons: Season identifiers\nharvests: Harvest identifiers\nsites: Site identifiers\nentries: Entry identifiers\npopulations: Population identifiers\nreplications: Replication identifiers\nblocks: Block identifiers\nrows: Row identifiers\ncols: Column identifiers\n\nAll remaining columns are treated as numeric phenotype measurements. Column names are fuzzy-matched to accommodate slight spelling variations.\n\nReturns\n\nTrials: A populated Trials struct containing the loaded data\n\nNotes\n\nMissing values can be represented as \"missing\", \"NA\", \"na\", \"N/A\", \"n/a\", or empty strings\nTrait names containing mathematical operators (+, -, *, /, %) are converted to underscores\nDuplicate trait names are not allowed\n\nThrows\n\nErrorException: If the input file doesn't exist or has invalid format\nArgumentError: If required columns are missing or ambiguous\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, sparsity=0.1, verbose=false);\n\njulia> fname = writedelimited(trials);\n\njulia> trials_reloaded = readdelimited(Trials, fname=fname);\n\njulia> trials == trials_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.readjld2-Union{Tuple{Type{T}}, Tuple{T}} where T<:GenomicBreedingCore.AbstractGB","page":"Home","title":"GenomicBreedingIO.readjld2","text":"readjld2(type::Type; fname::String)::Type\n\nLoad a core (Genomes, Phenomes, and Trials), simulation (SimulatedEffects), or model (TEBV) struct from a JLD2 file.\n\nArguments\n\ntype::Type: The type of struct to load (Genomes, Phenomes, Trials, SimulatedEffects, or TEBV)\nfname::String: Path to the JLD2 file to read from\n\nReturns\n\nThe loaded struct of the specified type\n\nThrows\n\nArgumentError: If the specified file does not exist\nDimensionMismatch: If the loaded struct is corrupted\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=2, verbose=false);\n\njulia> fname = writejld2(genomes);\n\njulia> readjld2(Genomes, fname=fname) == genomes\ntrue\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> fname = writejld2(phenomes);\n\njulia> readjld2(Phenomes, fname=fname) == phenomes\ntrue\n\njulia> trials, _ = simulatetrials(genomes=genomes, verbose=false);\n\njulia> fname = writejld2(trials);\n\njulia> readjld2(Trials, fname=fname) == trials\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.readvcf-Tuple{}","page":"Home","title":"GenomicBreedingIO.readvcf","text":"readvcf(; fname::String, field::String = \"any\", verbose::Bool = false)::Genomes\n\nRead genetic data from a VCF (Variant Call Format) file into a Genomes struct.\n\nArguments\n\nfname::String: Path to the VCF file. Can be gzipped (.vcf.gz or .vcf.bgz) or uncompressed (.vcf)\nfield::String=\"any\": Which FORMAT field to extract from VCF. Default \"any\" tries to automatically detect genotype field\nverbose::Bool=false: Whether to print progress and debug information\n\nReturns\n\nGenomes: A Genomes struct containing the loaded genetic data with fields:\nallele_frequencies: Matrix of allele frequencies\nloci_alleles: Vector of locus-allele combination strings\nmask: Boolean matrix indicating missing data\nsamples: Vector of sample names\n\nDetails\n\nReads VCF files in parallel using multiple threads. Handles multi-allelic variants and different ploidies.  Performs various checks on the input data including:\n\nFile existence\nNo duplicate loci-allele combinations \nConsistent dimensions in output struct\n\nThrows\n\nErrorException: If file doesn't exist, has duplicates, or output dimensions are invalid\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> fname_gz = writevcf(genomes, gzip=true);\n\njulia> genomes_reloaded = readvcf(fname=fname);\n\njulia> genomes_reloaded_gz = readvcf(fname=fname_gz);\n\njulia> genomes.entries == genomes_reloaded.entries == genomes_reloaded_gz.entries\ntrue\n\njulia> dimensions(genomes) == dimensions(genomes_reloaded) == dimensions(genomes_reloaded_gz)\ntrue\n\njulia> ismissing.(genomes.allele_frequencies) == ismissing.(genomes_reloaded.allele_frequencies) == ismissing.(genomes_reloaded_gz.allele_frequencies)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfchunkify-Tuple{String}","page":"Home","title":"GenomicBreedingIO.vcfchunkify","text":"vcfchunkify(fname::String; n_loci::Int64, verbose::Bool = false)::Tuple{Vector{Int64},Vector{Int64},Vector{Int64},Vector{Int64}}\n\nDivide a VCF file into chunks for parallel processing.\n\nArguments\n\nfname::String: Path to the VCF file (can be .vcf, .vcf.gz, or .vcf.bgz)\nn_loci::Int64: Total number of loci in the VCF file\nverbose::Bool=false: If true, prints progress information\n\nReturns\n\nA tuple containing four Vector{Int64} arrays:\n\nStarting loci indices for each thread\nEnding loci indices for each thread\nStarting file positions for each thread\nEnding file positions for each thread\n\nDetails\n\nAutomatically detects if the input file is gzipped\nDivides the workload evenly across available threads\nSkips header lines (starting with '#')\nHandles both regular and gzipped VCF files\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> _, n_loci = vcfcountloci(fname);\n\njulia> idx_loci_per_thread_ini, idx_loci_per_thread_fin, file_pos_per_thread_ini, file_pos_per_thread_fin = vcfchunkify(fname, n_loci=n_loci);\n\njulia> length(idx_loci_per_thread_ini) == length(idx_loci_per_thread_fin) == length(file_pos_per_thread_ini) == length(file_pos_per_thread_fin)\ntrue\n\njulia> (idx_loci_per_thread_ini[1] == 0) && (sum(idx_loci_per_thread_ini .== 0) == 1)\ntrue\n\njulia> (idx_loci_per_thread_fin[end] == n_loci) && (sum(idx_loci_per_thread_fin .== 0) == 0)\ntrue\n\njulia> (sum(file_pos_per_thread_ini .== 0) == 0) && (sum(file_pos_per_thread_fin .== 0) == 0)\ntrue\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfcountloci-Tuple{String}","page":"Home","title":"GenomicBreedingIO.vcfcountloci","text":"vcfcountloci(fname::String; verbose::Bool = false)::Tuple{Int64,Int64}\n\nCount the number of loci and total lines in a VCF file.\n\nArguments\n\nfname::String: Path to the VCF file. Can be either a plain text VCF file or a gzipped VCF file (with extensions .vcf.gz or .vcf.bgz)\nverbose::Bool: If true, prints progress messages and results to stdout. Defaults to false.\n\nReturns\n\nTuple{Int64,Int64}: A tuple containing:\nFirst element: Total number of lines in the file (including headers)\nSecond element: Number of data lines (variants/loci) excluding header lines\n\nDescription\n\nReads through a VCF (Variant Call Format) file and counts:\n\nTotal lines in the file (including headers)\nNumber of data lines (variants/loci) that don't start with '#'\n\nThe function automatically detects and handles different file formats:\n\nPlain text VCF files (.vcf)\nGzipped VCF files (.vcf.gz)\nBGZipped VCF files (.vcf.bgz)\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> fname_gz = writevcf(genomes, gzip=true);\n\njulia> n_1, p_1 = vcfcountloci(fname);\n\njulia> n_2, p_2 = vcfcountloci(fname_gz);\n\njulia> n_1 == n_2 == 10_009\ntrue\n\njulia> p_1 == p_2 == 10_000\ntrue\n\njulia> rm.([fname, fname_gz]);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfextractallelefreqs!-Tuple{GenomicBreedingCore.Genomes, Union{Nothing, ProgressMeter.Progress}, Vector{Int64}}","page":"Home","title":"GenomicBreedingIO.vcfextractallelefreqs!","text":"vcfextractallelefreqs!(input::Vector{Any}; fname::String, line::Vector{String}, \n                      line_counter::Int64, field::String, verbose::Bool = false)\n\nExtract allele frequencies from VCF file data and update a Genomes object.\n\nArguments\n\ninput::Vector{Any}: Vector containing:\ngenomes::Genomes: Object to store genomic data \npb::Union{Nothing,Progress}: Progress bar object or nothing\ni::Int64: Current locus-allele index\nfname::String: Name of VCF file being processed\nline::Vector{String}: Current line from VCF file split into fields\nline_counter::Int64: Current line number in VCF file\nfield::String: Type of field to extract (\"AF\", \"AD\", or \"GT\")\nverbose::Bool=false: Whether to display progress updates\n\nReturns\n\nUpdates the input vector's elements in place:\n\nFirst element (genomes): Updated with new allele frequencies\nSecond element (progress bar): Advanced if verbose=true\nThird element (index): Incremented based on processed alleles\n\nDescription\n\nProcesses VCF data to extract allele frequencies using:\n\nAF field: Direct frequency values\nAD field: Calculated from read depths\nGT field: Calculated from genotype calls \n\nUpdates Genomes object with coordinates and frequencies for each allele.\n\nThrows\n\nArgumentError: If input vector elements have incorrect types\nErrorException: If unable to parse AF or AD fields\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> _, n_loci = vcfcountloci(fname);\n\njulia> entries, format_lines = vcfextractentriesandformats(fname);\n\njulia> field, n_alleles, _ = vcfextractinfo(fname, format_lines=format_lines);\n\njulia> genomes_instantiated = vcfinstantiateoutput(fname, entries=entries, n_loci=n_loci, n_alleles=n_alleles);\n\njulia> sum(ismissing.(genomes_instantiated.allele_frequencies[:, 1])) == length(entries)\ntrue\n\njulia> file = open(fname, \"r\"); line::Vector{String} = split([readline(file) for i in 1:10][end], \"\t\"); close(file);\n\njulia> vcfextractallelefreqs!(genomes_instantiated, nothing, [0], fname=fname, line=line, line_counter=10, field=field);\n\njulia> sum(ismissing.(genomes_instantiated.allele_frequencies[:, 1])) == 0\ntrue\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfextractentriesandformats-Tuple{String}","page":"Home","title":"GenomicBreedingIO.vcfextractentriesandformats","text":"vcfextractentriesandformats(fname::String; verbose::Bool = false)::Tuple{Vector{String},Vector{String}}\n\nExtract sample entries and format definitions from a VCF file.\n\nArguments\n\nfname::String: Path to the VCF file (can be gzipped with extensions .vcf.gz or .vcf.bgz)\nverbose::Bool=false: If true, prints progress information to stdout\n\nReturns\n\nA tuple containing:\n\nVector{String}: List of sample names from the VCF header\nVector{String}: List of FORMAT field definitions from the VCF metadata\n\nDescription\n\nReads a VCF file and extracts two key pieces of information:\n\nSample names from the header line (columns after FORMAT field)\nFORMAT field definitions from metadata lines starting with \"##FORMAT\"\n\nThe function validates the presence and correct order of standard VCF columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT\n\nThrows\n\nArgumentError: If VCF has fewer than expected columns or column names don't match VCF format\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> entries, format_lines = vcfextractentriesandformats(fname);\n\njulia> entries\n10-element Vector{String}:\n \"entry_01\"\n \"entry_02\"\n \"entry_03\"\n \"entry_04\"\n \"entry_05\"\n \"entry_06\"\n \"entry_07\"\n \"entry_08\"\n \"entry_09\"\n \"entry_10\"\n\njulia> format_lines\n3-element Vector{String}:\n \"##FORMAT=<ID=GT,Number=1,Type=String,Description=\\\"Genotype\\\">\"\n \"##FORMAT=<ID=AD,Number=2,Type=Float,Description=\\\"Allele Depth\\\">\"\n \"##FORMAT=<ID=AF,Number=2,Type=Float,Description=\\\"Allele Frequency\\\">\"\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfextractinfo-Tuple{String}","page":"Home","title":"GenomicBreedingIO.vcfextractinfo","text":"vcfextractinfo(fname::String; format_lines::Vector{String}, field::String=\"any\", verbose::Bool=false)::Tuple{String,Int64,Int64}\n\nExtract information about genotype fields from a VCF file.\n\nArguments\n\nfname::String: Path to the VCF file (can be gzipped)\nformat_lines::Vector{String}: Vector containing FORMAT lines from the VCF header\nfield::String=\"any\": Specific field to extract (\"GT\", \"AD\", \"AF\", or \"any\")\nverbose::Bool=false: If true, prints progress information\n\nReturns\n\nA tuple containing:\n\nfield::String: The identified genotype field\nn_alleles::Int64: Maximum number of alleles per locus\nploidy::Int64: Ploidy level (only meaningful for GT field; set to typemax(Int64) for AD and AF fields)\n\nDetails\n\nIf field is \"any\", searches for fields in priority order: AF > AD > GT\nFor GT field, scans entire file to determine maximum number of alleles and ploidy\nFor AF and AD fields, extracts allele count from format header\nSupports both gzipped (.gz, .bgz) and uncompressed VCF files\n\nThrows\n\nArgumentError: If specified field is not found in the VCF file\nErrorException: If unable to parse number of alleles from format header\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> _, format_lines = vcfextractentriesandformats(fname);\n\njulia> field, n_alleles, ploidy = vcfextractinfo(fname, format_lines=format_lines);\n\njulia> (field == \"AF\") && (n_alleles == 2) && (ploidy == typemax(Int64))\ntrue\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfinstantiateoutput-Tuple{String}","page":"Home","title":"GenomicBreedingIO.vcfinstantiateoutput","text":"vcfinstantiateoutput(fname::String; entries::Vector{String}, n_loci::Int64, n_alleles::Int64, verbose::Bool = false)::Genomes\n\nCreate and initialize a Genomes struct from VCF file parameters.\n\nArguments\n\nfname::String: Name of the VCF file being processed\nentries::Vector{String}: Vector containing entry identifiers\nn_loci::Int64: Number of loci in the VCF file\nn_alleles::Int64: Number of alleles per locus\nverbose::Bool=false: If true, prints progress information\n\nReturns\n\nGenomes: An initialized Genomes struct with:\ndimensions n × p where n is number of entries and p = nloci * (nalleles - 1)\nentry names assigned\npopulations set to \"unknown\"\nmask set to true\n\nThrows\n\nErrorException: If duplicate entries are found in the VCF file\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> entries, format_lines = vcfextractentriesandformats(fname);\n\njulia> _, n_loci = vcfcountloci(fname);\n\njulia> _, n_alleles, _ = vcfextractinfo(fname, format_lines=format_lines);\n\njulia> genomes_instantiated = vcfinstantiateoutput(fname, entries=entries, n_loci=n_loci, n_alleles=n_alleles);\n\njulia> size(genomes_instantiated.allele_frequencies)\n(10, 10000)\n\njulia> genomes_instantiated.entries == entries\ntrue\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.vcfparsecoordinates-Tuple{}","page":"Home","title":"GenomicBreedingIO.vcfparsecoordinates","text":"vcfparsecoordinates(; line::Vector{String}, line_counter::Int64, field::String)::Union{Nothing,Tuple{Int64,String,Int64,Vector{String}}}\n\nParse coordinates and allele information from a VCF file line.\n\nArguments\n\nline::Vector{String}: A vector containing the split line from VCF file\nline_counter::Int64: Current line number being processed in the VCF file\nfield::String: The field name to extract allele frequencies from\n\nReturns\n\nNothing: If the specified field is not found in the line\nTuple{Int64,String,Int64,Vector{String}}: A tuple containing:\nField index\nChromosome name\nPosition\nCombined reference and alternative alleles\n\nThrows\n\nErrorException: If the position field cannot be parsed as an integer\n\nNote\n\nThe function validates the line format and extracts genomic coordinates and allele information from a VCF file line. It handles missing alternative alleles (denoted by \".\") and performs necessary type conversions.\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> entries, format_lines = vcfextractentriesandformats(fname);\n\njulia> field, _, _ = vcfextractinfo(fname, format_lines=format_lines);\n\njulia> file = open(fname, \"r\"); line::Vector{String} = split([readline(file) for i in 1:10][end], \"\t\"); close(file);\n\njulia> idx_field, chrom, pos, refalts = vcfparsecoordinates(line=line, line_counter=10, field=field);\n\njulia> (idx_field == 3) && (chrom == line[1]) && (pos == parse(Int64, line[2])) && (refalts == line[4:5])\ntrue\n\njulia> rm(fname);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.writedelimited-Tuple{GenomicBreedingCore.Genomes}","page":"Home","title":"GenomicBreedingIO.writedelimited","text":"writedelimited(\n    genomes::Genomes;\n    fname::Union{Missing,String} = missing,\n    sep::String = \"\\t\",\n    include_population_header::Bool = true\n)::String\n\nWrite genomic data to a delimited text file.\n\nArguments\n\ngenomes::Genomes: A Genomes struct containing the genomic data to be written\nfname::Union{Missing,String}: Output filename. If missing, generates an automatic filename with timestamp\nsep::String: Delimiter character for the output file (default: tab)\ninclude_population_header::Bool: Whether to include population information in the header (default: true)\n\nReturns\n\nString: Path to the created output file\n\nFile Format\n\nThe output file contains:\n\nHeader lines (prefixed with '#'):\nFirst line: chromosome, position, alleles, and entry information\nSecond line (optional): population information\nData rows with the following columns:\nColumn 1: Chromosome identifier\nColumn 2: Position\nColumn 3: All alleles at the locus (pipe-separated)\nColumn 4: Specific allele\nRemaining columns: Frequency data for each entry\n\nSupported File Extensions\n\n'.tsv' (tab-separated, default)\n'.csv' (comma-separated)\n'.txt' (custom delimiter)\n\nThrows\n\nDimensionMismatch: If the input Genomes struct is corrupted\nErrorException: If the output file already exists\nArgumentError: If the file extension is invalid or the output directory doesn't exist\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=2, verbose=false);\n\njulia> writedelimited(genomes, fname=\"test_genomes.tsv\")\n\"test_genomes.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.writedelimited-Tuple{GenomicBreedingCore.Phenomes}","page":"Home","title":"GenomicBreedingIO.writedelimited","text":"writedelimited(phenomes::Phenomes; fname::Union{Missing,String} = missing, sep::String = \"\t\")::String\n\nWrite phenotypic data from a Phenomes struct to a delimited text file.\n\nArguments\n\nphenomes::Phenomes: A Phenomes struct containing phenotypic data\nfname::Union{Missing,String} = missing: Output filename. If missing, generates an automatic filename with timestamp\nsep::String = \"\t\": Delimiter character for the output file\n\nReturns\n\nString: The name of the created file\n\nFile Format\n\nHeader line starts with '#' containing column names\nFirst column: Entry names\nSecond column: Population names\nRemaining columns: Trait values\nMissing values are represented as \"NA\"\n\nFile Extensions\n\nSupported file extensions:\n\n.tsv for tab-separated files (default)\n.csv for comma-separated files\n.txt for other delimiters\n\nThrows\n\nDimensionMismatch: If the Phenomes struct dimensions are inconsistent\nErrorException: If the output file already exists\nArgumentError: If the file extension is invalid or the directory doesn't exist\n\nExamples\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writedelimited(phenomes, fname=\"test_phenomes.tsv\")\n\"test_phenomes.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.writedelimited-Tuple{GenomicBreedingCore.Trials}","page":"Home","title":"GenomicBreedingIO.writedelimited","text":"writedelimited(trials::Trials; fname::Union{Missing,String} = missing, sep::String = \"\t\")::String\n\nWrite a Trials struct to a delimited text file, returning the filename.\n\nArguments\n\ntrials::Trials: The trials data structure to be written\nfname::Union{Missing,String} = missing: Output filename. If missing, generates automatic filename with timestamp\nsep::String = \"\t\": Delimiter character between fields\n\nReturns\n\nString: The name of the file that was written\n\nFile Format\n\nThe output file contains one header line and one line per trial entry. Header line is prefixed with '#' and contains column names.\n\nFixed Columns (1-10)\n\nyears\nseasons\nharvests\nsites\nentries\npopulations\nreplications\nblocks\nrows\ncols\n\nVariable Columns (11+)\n\nAdditional columns contain phenotype traits values\nMissing values are written as \"NA\"\n\nNotes\n\nSupported file extensions: .tsv, .csv, or .txt\nFile extension is automatically determined based on separator if filename is missing:\n\\t → .tsv\n, or ; → .csv\nother → .txt\nWill not overwrite existing files\nDirectory must exist if path is specified in filename\n\nExamples\n\njulia> trials = Trials(n=1, t=2); trials.years = [\"year_1\"]; trials.seasons = [\"season_1\"]; trials.harvests = [\"harvest_1\"]; trials.sites = [\"site_1\"]; trials.entries = [\"entry_1\"]; trials.populations = [\"population_1\"]; trials.replications = [\"replication_1\"]; trials.blocks = [\"block_1\"]; trials.rows = [\"row_1\"]; trials.cols = [\"col_1\"]; trials.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writedelimited(trials, fname=\"test_trials.tsv\")\n\"test_trials.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.writejld2-Tuple{GenomicBreedingCore.AbstractGB}","page":"Home","title":"GenomicBreedingIO.writejld2","text":"writejld2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects,TEBV}; fname::Union{Missing,String} = missing)::String\n\nSave genomic breeding core data structures to a JLD2 file (HDF5-compatible format).\n\nArguments\n\nA: A genomic breeding data structure (Genomes, Phenomes, Trials, SimulatedEffects, or TEBV)\nfname: Optional. Output filename. If missing, generates an automatic name with timestamp\n\nReturns\n\nString: Path to the saved JLD2 file\n\nFile Naming\n\nIf fname is not provided, generates name: \"output-[Type]-[Timestamp].jld2\"\nIf fname is provided, must have \".jld2\" extension\n\nThrows\n\nDimensionMismatch: If input structure has invalid dimensions\nErrorException: If output file already exists\nArgumentError: If invalid file extension or directory path\n\nNotes\n\nFiles are saved with compression enabled\nData is stored as a Dictionary with single key-value pair\nKey is the string representation of the input type\nExisting files will not be overwritten\n\nExamples\n\njulia> genomes = GenomicBreedingCore.simulategenomes(n=2, verbose=false);\n\njulia> writejld2(genomes, fname=\"test_genomes.jld2\")\n\"test_genomes.jld2\"\n\njulia> genomes_reloaded = load(\"test_genomes.jld2\");\n\njulia> genomes_reloaded[collect(keys(genomes_reloaded))[1]] == genomes\ntrue\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writejld2(phenomes, fname=\"test_phenomes.jld2\")\n\"test_phenomes.jld2\"\n\njulia> phenomes_reloaded = load(\"test_phenomes.jld2\");\n\njulia> phenomes_reloaded[collect(keys(phenomes_reloaded))[1]] == phenomes\ntrue\n\njulia> trials, _ = simulatetrials(genomes=genomes, verbose=false);\n\njulia> writejld2(trials, fname=\"test_trials.jld2\")\n\"test_trials.jld2\"\n\njulia> trials_reloaded = load(\"test_trials.jld2\");\n\njulia> trials_reloaded[collect(keys(trials_reloaded))[1]] == trials\ntrue\n\njulia> simulated_effects = SimulatedEffects();\n\njulia> writejld2(simulated_effects, fname=\"test_simulated_effects.jld2\")\n\"test_simulated_effects.jld2\"\n\njulia> simulated_effects_reloaded = load(\"test_simulated_effects.jld2\");\n\njulia> simulated_effects_reloaded[collect(keys(simulated_effects_reloaded))[1]] == simulated_effects\ntrue\n\njulia> trials, _simulated_effects = GenomicBreedingCore.simulatetrials(genomes = GenomicBreedingCore.simulategenomes(n=10, verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false);\n\njulia> tebv = analyse(trials, max_levels=50, verbose=false);\n\njulia> writejld2(tebv, fname=\"test_tebv.jld2\")\n\"test_tebv.jld2\"\n\njulia> tebv_reloaded = load(\"test_tebv.jld2\");\n\njulia> tebv_reloaded[collect(keys(tebv_reloaded))[1]] == tebv\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreedingIO.writevcf-Tuple{GenomicBreedingCore.Genomes}","page":"Home","title":"GenomicBreedingIO.writevcf","text":"writevcf(genomes::Genomes; fname::Union{Missing,String} = missing, ploidy::Int64 = 0, \n         max_depth::Int64 = 100, n_decimal_places::Int64 = 4, gzip::Bool = false)::String\n\nWrite genomic data to a Variant Call Format (VCF) file.\n\nArguments\n\ngenomes::Genomes: A Genomes object containing the genetic data to be written.\nfname::Union{Missing,String} = missing: Output filename. If missing, generates a default name with timestamp.\nploidy::Int64 = 0: The ploidy level of the organisms (e.g., 2 for diploid).\nmax_depth::Int64 = 100: Maximum depth for allele depth calculation.\nn_decimal_places::Int64 = 4: Number of decimal places for rounding allele frequencies.\ngzip::Bool = false: Whether to compress the output file using gzip.\n\nReturns\n\nString: The name of the created VCF file.\n\nDescription\n\nCreates a VCF v4.2 format file containing genomic variants data. The function processes allele frequencies and depths, calculates genotypes based on ploidy, and formats the data according to VCF specifications. The output includes:\n\nStandard VCF header information\nSample information with FORMAT fields:\nGT (Genotype)\nAD (Allele Depth)\nAF (Allele Frequency)\n\nThrows\n\nDimensionMismatch: If the input Genomes object has inconsistent dimensions\nErrorException: If the output file already exists\nArgumentError: If the file extension is not '.vcf' or if the specified directory doesn't exist\n\nExamples\n\njulia> genomes_1 = GenomicBreedingCore.simulategenomes(n=2, verbose=false);\n\njulia> writevcf(genomes_1, fname=\"test_genomes_1.vcf\")\n\"test_genomes_1.vcf\"\n\njulia> genomes_2 = GenomicBreedingCore.simulategenomes(n=2, n_alleles=3, verbose=false);\n\njulia> genomes_2.allele_frequencies = round.(genomes_2.allele_frequencies .* 4) ./ 4;\n\njulia> writevcf(genomes_2, fname=\"test_genomes_2.vcf\", ploidy=4)\n\"test_genomes_2.vcf\"\n\njulia> genomes_3 = GenomicBreedingCore.simulategenomes(n=3, verbose=false);\n\njulia> writevcf(genomes_3, fname=\"test_genomes_3.vcf\", gzip=true)\n\"test_genomes_3.vcf.gz\"\n\n\n\n\n\n","category":"method"}]
}
