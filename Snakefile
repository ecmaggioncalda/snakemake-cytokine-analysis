configfile: 'config/config.yml'

core_genome_data = config['core_genome_data']
pan_genome_data = config['pan_genome_data']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']
outcome_colname = config['outcome_colname']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule prepro_overall:
    input:
        R = "code/prepro_overall.R",
        patient_metadata = config['patient_metadata']
    output:
        "data/{cytokine}.tsv"
        "data/{cytokine}_adjusted.tsv"
    log: "log/prepro_overall.txt"
    params:
    resources:
        ncores = ncores
    scripts:
        "code/prepro_overall.R"

checkpoint adjustments:
    input:
        R = "code/verify_adjustment_files.R"
        "data/{cytokine}.tsv"
    output:
        directory("/{cytokine}_adjusted")
    scripts:
        "code/verify_adjustment_files.R"

rule preprocess_frames:
    input:
    output:
    log:
    params:
    resources:
    scripts:

rule core_genome_hogwash:
    input: "combined_core.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule core_genome_elastic_net:
    input: "combined_core.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule core_genome_treeWAS:
    input: "combined_core.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule pan_genome_hogwash:
    input: "pan_matrix.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule pan_genome_elastic_net:
    input: "pan_matrix.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule pan_genome_treeWAS:
    input: "pan_matrix.tsv"
    output:
    log:
    params:
    resources:
    scripts:

rule summarize_cytokine_results:
    input:
    output:
    log:
    params:
    resources:
    scripts:
