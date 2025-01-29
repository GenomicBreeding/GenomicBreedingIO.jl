var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GBIO","category":"page"},{"location":"#GBIO","page":"Home","title":"GBIO","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GBIO.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GBIO]","category":"page"},{"location":"#GBIO.isfuzzymatch-Tuple{String, String}","page":"Home","title":"GBIO.isfuzzymatch","text":"isfuzzymatch(a::String, b::String; threshold::Float64=0.3)::Bool\n\nFuzzy string matching using the Levenshtein distance, and a threshold as a fraction of the smaller string.\n\nExamples\n\njulia> isfuzzymatch(\"populations\", \"populations\")\ntrue\n\njulia> isfuzzymatch(\"populations\", \"poplation\")\ntrue\n\njulia> isfuzzymatch(\"populations\", \"entry\")\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.levenshteindistance-Tuple{String, String}","page":"Home","title":"GBIO.levenshteindistance","text":"levenshteindistance(a::String, b::String)::Int64\n\nCalculate the Levenshtein distance between 2 strings. TODO: optimise as we only need to compute 2 rows at a time and do not need the full matrix\n\nExamples\n\njulia> levenshteindistance(\"populations\", \"populations\")\n0\n\njulia> levenshteindistance(\"populations\", \"poplation\")\n2\n\njulia> levenshteindistance(\"populations\", \"entry\")\n3\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.readJLD2-Union{Tuple{Type{T}}, Tuple{T}, Tuple{Type{T}, String}} where T<:GBCore.AbstractGB","page":"Home","title":"GBIO.readJLD2","text":"readJLD2(type::Type, fname::String = missing)::Type\n\nLoad a core (Genomes, Phenomes, and Trials), simulation (SimulatedEffects), and model (TEBV) struct from a JLD2 file.\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=2, verbose=false);\n\njulia> fname = writeJLD2(genomes);\n\njulia> readJLD2(Genomes, fname) == genomes\ntrue\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> fname = writeJLD2(phenomes);\n\njulia> readJLD2(Phenomes, fname) == phenomes\ntrue\n\njulia> trials, _ = simulatetrials(genomes=genomes, verbose=false);\n\njulia> fname = writeJLD2(trials);\n\njulia> readJLD2(Trials, fname) == trials\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.readdelimited-Tuple{Type{GBCore.Genomes}}","page":"Home","title":"GBIO.readdelimited","text":"readdelimited(type::Type{Genomes}; fname::String, sep::String = \"\\t\")::Genomes\n\nLoad a Genomes struct from a string-delimited (default=\\t) file.  Each row corresponds to a locus-allele combination. The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by |), and the specific allele. The subsequency columns refer to the samples, pools, entries or genotypes.\n\nNotes:\n\nExtension name should be '.tsv', '.csv', or '.txt'.\nHeader lines and comments are prefixed by '#'.\nThere are 2 header lines prefixed by '#', e.g.:\nheader line 1: \"chrom,pos,allalleles,allele,entry1,entry_2\"\nheader line 2: \"chrom,pos,allalleles,allele,population1,population_1\"\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writedelimited(genomes);\n\njulia> genomes_reloaded = readdelimited(Genomes, fname=fname);\n\njulia> genomes == genomes_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.readdelimited-Tuple{Type{GBCore.Phenomes}}","page":"Home","title":"GBIO.readdelimited","text":"readdelimited(type::Type{Phenomes}; fname::String, sep::String = \"\\t\")::Phenomes\n\nLoad a Phenomes struct from a string-delimited (default=\\t) file.  Each row corresponds to a locus-allele combination. The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by |), and the specific allele. The subsequency columns refer to the samples, pools, entries or genotypes.\n\nExamples\n\njulia> phenomes = Phenomes(n=10, t=3); phenomes.entries = string.(\"entry_\", 1:10); phenomes.populations .= \"pop1\"; phenomes.traits = [\"A\", \"B\", \"C\"]; phenomes.phenotypes = rand(10,3); phenomes.phenotypes[1,1] = missing; phenomes.mask .= true;\n\njulia> fname = writedelimited(phenomes);\n\njulia> phenomes_reloaded = readdelimited(Phenomes, fname=fname);\n\njulia> phenomes == phenomes_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.readdelimited-Tuple{Type{GBCore.Trials}}","page":"Home","title":"GBIO.readdelimited","text":"readdelimited(type::Type{Trials}; fname::String, sep::String = \"\\t\")::Trials\n\nLoad a Trials struct from a string-delimited (default=\\t) file.  We expect the following 10 identifier columns:     ‣ years     ‣ seasons     ‣ harvests     ‣ sites     ‣ entries     ‣ populations     ‣ replications     ‣ blocks     ‣ rows     ‣ cols All other columns are assumed to be numeric phenotype values.\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=10, verbose=false);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, sparsity=0.1, verbose=false);\n\njulia> fname = writedelimited(trials);\n\njulia> trials_reloaded = readdelimited(Trials, fname=fname);\n\njulia> trials == trials_reloaded\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.readvcf-Tuple{}","page":"Home","title":"GBIO.readvcf","text":"readvcf(;fname::String, field::String = \"any\", verbose::Bool = false)::Genomes\n\nLoad Genomes struct from vcf file\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=10, verbose=false);\n\njulia> fname = writevcf(genomes);\n\njulia> fname_gz = writevcf(genomes, gzip=true);\n\njulia> genomes_reloaded = readvcf(fname=fname);\n\njulia> genomes_reloaded_gz = readvcf(fname=fname_gz);\n\njulia> genomes.entries == genomes_reloaded.entries == genomes_reloaded_gz.entries\ntrue\n\njulia> dimensions(genomes) == dimensions(genomes_reloaded) == dimensions(genomes_reloaded_gz)\ntrue\n\njulia> ismissing.(genomes.allele_frequencies) == ismissing.(genomes_reloaded.allele_frequencies) == ismissing.(genomes_reloaded_gz.allele_frequencies)\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.writeJLD2-Tuple{GBCore.AbstractGB}","page":"Home","title":"GBIO.writeJLD2","text":"writeJLD2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects}; fname::Union{Missing,String} = missing)::String\n\nSave core (Genomes, Phenomes, and Trials), simulation (SimulatedEffects), and model (TEBV) structs as JLD2,  a heirarchical data format version 5 (HDF5) - compatible format. Note that the extension name should be '.jld2'.\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=2, verbose=false);\n\njulia> writeJLD2(genomes, fname=\"test_genomes.jld2\")\n\"test_genomes.jld2\"\n\njulia> genomes_reloaded = load(\"test_genomes.jld2\");\n\njulia> genomes_reloaded[collect(keys(genomes_reloaded))[1]] == genomes\ntrue\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writeJLD2(phenomes, fname=\"test_phenomes.jld2\")\n\"test_phenomes.jld2\"\n\njulia> phenomes_reloaded = load(\"test_phenomes.jld2\");\n\njulia> phenomes_reloaded[collect(keys(phenomes_reloaded))[1]] == phenomes\ntrue\n\njulia> trials, _ = simulatetrials(genomes=genomes, verbose=false);\n\njulia> writeJLD2(trials, fname=\"test_trials.jld2\")\n\"test_trials.jld2\"\n\njulia> trials_reloaded = load(\"test_trials.jld2\");\n\njulia> trials_reloaded[collect(keys(trials_reloaded))[1]] == trials\ntrue\n\njulia> simulated_effects = SimulatedEffects();\n\njulia> writeJLD2(simulated_effects, fname=\"test_simulated_effects.jld2\")\n\"test_simulated_effects.jld2\"\n\njulia> simulated_effects_reloaded = load(\"test_simulated_effects.jld2\");\n\njulia> simulated_effects_reloaded[collect(keys(simulated_effects_reloaded))[1]] == simulated_effects\ntrue\n\njulia> trials, _simulated_effects = GBCore.simulatetrials(genomes = GBCore.simulategenomes(n=10, verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false);\n\njulia> tebv = analyse(trials, max_levels=50, verbose=false);\n\njulia> writeJLD2(tebv, fname=\"test_tebv.jld2\")\n\"test_tebv.jld2\"\n\njulia> tebv_reloaded = load(\"test_tebv.jld2\");\n\njulia> tebv_reloaded[collect(keys(tebv_reloaded))[1]] == tebv\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.writedelimited-Tuple{GBCore.Genomes}","page":"Home","title":"GBIO.writedelimited","text":"writedelimited(genomes::Genomes, sep::String = \"\t\", fname::Union{Missing,String} = missing)::String\n\nSave Genomes struct as a string-delimited (default=) file. Each row corresponds to a locus-allele combination. The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by |), and the specific allele. The subsequency columns refer to the samples, pools, entries or genotypes.\n\nNotes:\n\nExtension name should be '.tsv', '.csv', or '.txt'.\nHeader lines and comments are prefixed by '#'.\nThere are 2 header lines prefixed by '#', e.g.:\nheader line 1: \"chrom,pos,allalleles,allele,entry1,entry_2\"\nheader line 2: \"chrom,pos,allalleles,allele,population1,population_1\"\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=2, verbose=false);\n\njulia> writedelimited(genomes, fname=\"test_genomes.tsv\")\n\"test_genomes.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.writedelimited-Tuple{GBCore.Phenomes}","page":"Home","title":"GBIO.writedelimited","text":"writedelimited(phenomes::Phenomes, sep::String = \"\t\", fname::Union{Missing,String} = missing)::String\n\nSave Phenomes struct as a string-delimited (default=) file.  Each row corresponds to a samples, pools, entries or genotypes. The first 2 columns correspond to the entry and population names. The subsequency columns refer to the traits containing the phenotype values of each entry. Note that the extension name should be '.tsv', '.csv', or '.txt'.\n\nNotes:\n\nExtension name should be '.tsv', '.csv', or '.txt'.\nHeader line and comments are prefixed by '#'.\nThere is 1 header line prefixed by '#', e.g.: \"entry,population,trait1,trait2,trait_3\"\n\nExamples\n\njulia> phenomes = Phenomes(n=2, t=2); phenomes.entries = [\"entry_1\", \"entry_2\"]; phenomes.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writedelimited(phenomes, fname=\"test_phenomes.tsv\")\n\"test_phenomes.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.writedelimited-Tuple{GBCore.Trials}","page":"Home","title":"GBIO.writedelimited","text":"writedelimited(trials::Trials, sep::String = \"\t\", fname::Union{Missing,String} = missing)::String\n\nSave Trials struct as a string-delimited (default=) file.  Each row corresponds to a samples, pools, entries or genotypes. The first 10 columns correspond to:\n\nyears\nseasons\nharvests\nsites\nentries\npopulations\nreplications\nblocks\nrows\ncols \n\nThe subsequency columns refer to the traits containing the phenotype values.\n\nNotes:\n\nExtension name should be '.tsv', '.csv', or '.txt'.\nHeader line and comments are prefixed by '#'.\nThere is 1 header line prefixed by '#', e.g.: \"years,seasons,harvests, ..., trait1,tratit2,trait_3\"\n\nExamples\n\njulia> trials = Trials(n=1, t=2); trials.years = [\"year_1\"]; trials.seasons = [\"season_1\"]; trials.harvests = [\"harvest_1\"]; trials.sites = [\"site_1\"]; trials.entries = [\"entry_1\"]; trials.populations = [\"population_1\"]; trials.replications = [\"replication_1\"]; trials.blocks = [\"block_1\"]; trials.rows = [\"row_1\"]; trials.cols = [\"col_1\"]; trials.traits = [\"trait_1\", \"trait_2\"];\n\njulia> writedelimited(trials, fname=\"test_trials.tsv\")\n\"test_trials.tsv\"\n\n\n\n\n\n","category":"method"},{"location":"#GBIO.writevcf-Tuple{GBCore.Genomes}","page":"Home","title":"GBIO.writevcf","text":"writevcf(genomes::Genomes; fname::Union{Missing,String}, ploidy::Int64=0, max_depth::Int64=100, n_decimal_places::Int64=4)::String\n\nSave Genomes struct as a variant call format (VCF version 4.2) file.\n\nExamples\n\njulia> genomes_1 = GBCore.simulategenomes(n=2, verbose=false);\n\njulia> writevcf(genomes_1, fname=\"test_genomes_1.vcf\")\n\"test_genomes_1.vcf\"\n\njulia> genomes_2 = GBCore.simulategenomes(n=2, n_alleles=3, verbose=false);\n\njulia> genomes_2.allele_frequencies = round.(genomes_2.allele_frequencies .* 4) ./ 4;\n\njulia> writevcf(genomes_2, fname=\"test_genomes_2.vcf\", ploidy=4)\n\"test_genomes_2.vcf\"\n\njulia> genomes_3 = GBCore.simulategenomes(n=3, verbose=false);\n\njulia> writevcf(genomes_3, fname=\"test_genomes_3.vcf\", gzip=true)\n\"test_genomes_3.vcf.gz\"\n\n\n\n\n\n","category":"method"}]
}
