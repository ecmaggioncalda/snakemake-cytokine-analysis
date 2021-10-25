configfile: 'config/config.yml'

ncores = config['ncores']
cytokine = config['cytokine']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule determine_adjustments:
    input: "patient_metadata.csv"
    output: "{cytokine}_adjusted.tsv" "{cytokine}_unadjusted.tsv"
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
