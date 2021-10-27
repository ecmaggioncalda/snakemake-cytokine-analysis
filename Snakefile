configfile: 'config/config.yml'

tree = config['tree']
core = config['core']
pan = config['pan']
group = ["raw", "adjusted"]

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule prepro_overall:
    input:
        R = "code/prepro_overall.R",
        patient_metadata = config['patient_metadata']
    output:
        "data/{cytokine}_{group}.tsv"
    log:
        "log/prepro_overall.txt"
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
        R = "code/run_treeWAS.R"
        geno = config['pan']
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        rdata = 'results/{cytokine}/results/{cytokine}_{group}_pan_treeWAS.RData'
        plot = 'results/{cytokine}/plots/{cytokine}_{group}_pan_treeWAS.pdf'
    log:
        "log/{cytokine}_{group}_pan_treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped_pan:
    input:
        R = "code/run_hogwash_ungrouped.R"
        geno = config['pan']
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        file_name = '{cytokine}_{group}_pan_hogwash_ungrouped'
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}_{group}_pan_hogwash_ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

rule run_treeWAS_core:
    input:
        R = "code/run_treeWAS.R"
        geno = config['core']
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        rdata = 'results/{cytokine}/results/{cytokine}_{group}_core_treeWAS.RData'
        plot = 'results/{cytokine}/plots/{cytokine}_{group}_core_treeWAS.pdf'
    log:
        "log/{cytokine}_{group}_core_treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped_core:
    input:
        R = "code/run_hogwash_ungrouped.R"
        geno = config['core']
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        file_name = '{cytokine}_{group}_core_hogwash_ungrouped'
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}_{group}_core_hogwash_ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

rule run_hogwash_grouped:
    input:
        R = "code/run_hogwash_grouped.R"
        geno = config['core']
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
        gene_key = config['gene_key']
    output:
        file_name = '{cytokine}_{group}_core_hogwash_grouped'
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}_{group}_core_hogwash_grouped.txt"
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
