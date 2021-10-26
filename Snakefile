configfile: 'config/config.yml'

patient_metadata = config['patient_metadata']
core_genome_data = config['core_genome_data']
pan_genome_data = config['pan_genome_data']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']
outcome_colname = config['outcome_colname']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule determine_adjustments:
    input: patient_metadata
    output: "{cytokine}_adjusted.tsv" "{cytokine}_unadjusted.tsv"
    log:
    params:
    resources:
    scripts:

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
