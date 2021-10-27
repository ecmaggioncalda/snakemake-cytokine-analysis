configfile: 'config/config.yml'

tree = config['tree']
core = config['core']
pan = config['pan']
cytokine = config['cytokine']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

group = ["raw", "adjusted"]
type = ["pan", "core"]
model = ["hogwash_ungrouped", "treeWAS.RData"]

rule all_test:
    input:
        expand("results/{cytokine}/results/{group}_core_hogwash_grouped", cytokine = cytokine, group = group),
        expand("results/{cytokine}/results/{group}_{type}_{model}", cytokine = cytokine, group = group, type = type, model = model)

rule prepro_overall:
    input:
        R = "code/prepro_overall.R",
        patient_metadata = config['patient_metadata']
    output:
        "{cytokine}_{group}.tsv"
    params:
        cytokine = '{cytokine}'
    log:
        "log/{cytokine}/{group}_prepro_overall.txt"
    resources:
        ncores = ncores
    script:
        "code/prepro_overall.R"

###Come back to this, need a checkpoint to trim the initial DAG since not all cytokines have an adjusted file
#checkpoint adjustments:
#    input:
#        R = "code/verify_adjustment_files.R"
#        "data/{cytokine}.tsv"
#    output:
#
#    script:
#        "code/verify_adjustment_files.R"

rule run_treeWAS_pan:
    input:
        R = "code/run_treeWAS.R",
        pheno = 'data/{cytokine}_{group}.tsv'
    output:
        rdata = 'results/{cytokine}/results/{group}_pan_treeWAS.RData',
        plot = 'results/{cytokine}/plots/{group}_pan_treeWAS.pdf'
    params:
        geno = pan,
        tree = tree
    log:
        "log/{cytokine}/{group}_pan_treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped_pan:
    input:
        R = "code/run_hogwash_ungrouped.R",
        pheno = 'data/{cytokine}_{group}.tsv'
    output:
        "results/{cytokine}/results/{group}_pan_hogwash_ungrouped"
    params:
        file_name = '{cytokine}_pan_hogwash_ungrouped',
        geno = pan,
        tree = tree,
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}/{group}_pan_hogwash_ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

rule run_treeWAS_core:
    input:
        R = "code/run_treeWAS.R",
        pheno = 'data/{cytokine}_{group}.tsv'
    output:
        rdata = 'results/{cytokine}/results/{group}_core_treeWAS.RData',
        plot = 'results/{cytokine}/plots/{group}_core_treeWAS.pdf'
    params:
        geno = core,
        tree = tree
    log:
        "log/{cytokine}/{group}_core_treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped_core:
    input:
        R = "code/run_hogwash_ungrouped.R",
        pheno = 'data/{cytokine}_{group}.tsv'
    output:
        "results/{cytokine}/results/{group}_core_hogwash_ungrouped"
    params:
        file_name = '{group}_core_hogwash_ungrouped',
        geno = core,
        tree = tree,
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}/{group}_core_hogwash_ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

rule run_hogwash_grouped:
    input:
        R = "code/run_hogwash_grouped.R",
        pheno = 'data/{cytokine}_{group}.tsv'
    output:
        "results/{cytokine}/results/{group}_core_hogwash_grouped"
    params:
        file_name = '{cytokine}_core_hogwash_grouped',
        geno = core,
        tree = tree,
        dir = "results/{cytokine}/results/",
        gene_key = config['gene_key']
    log:
        "log/{cytokine}/{group}_core_hogwash_grouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_grouped.R"

#rule elastic_net:
#    input: "combined_core.tsv"
#    output:
#    log:
#    params:
#    resources:
#    script:

#rule summarize_cytokine_results:
#    input:
#    output:
#    log:
#    params:
#    resources:
#    script:
