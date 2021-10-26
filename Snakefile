configfile: 'config/config.yml'

tree = config['tree']
core = config['core_genome_data']
pan = config['pan_genome_data']
genome_data = [core, pan]
group = ["raw", "adjusted"]

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
        "data/{cytokine}_{group}.tsv"
        directory("results/{cytokine}")
        directory("results/{cytokine}/plots")
        directory("results/{cytokine}/results")
    log: "log/prepro_overall.txt"
    params:
    resources:
        ncores = ncores
    scripts:
        "code/prepro_overall.R"

###Come back to this, need a checkpoint to trim the initial DAG since not all cytokines have an adjusted file
#checkpoint adjustments:
#    input:
#        R = "code/verify_adjustment_files.R"
#        "data/{cytokine}.tsv"
#    output:
#
#    scripts:
#        "code/verify_adjustment_files.R"

rule run_treeWAS:
    input:
        R = "code/run_treeWAS.R"
        geno = genome_data
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        rdata = 'results/{cytokine}/results/{cytokine}_{group}_{geno}_treeWAS.RData'
        plot = 'results/{cytokine}/plots/{cytokine}_{group}_{geno}_treeWAS.pdf'
    log:
        "log/{cytokine}_{group}_{geno}_treeWAS.txt"
    resources:
        ncores = ncores
    scripts:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped:
    input:
        R = "code/run_hogwash_ungrouped.R"
        geno = genome_data
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
    output:
        file_name = '{cytokine}_{group}_{geno}_hogwash_ungrouped'
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}_{group}_{geno}_hogwash_ungrouped.txt"
    resources:
        ncores = ncores
    scripts:
        "code/run_hogwash_ungrouped.R"

rule run_hogwash_grouped:
    input:
        R = "code/run_hogwash_grouped.R"
        geno = core
        pheno = 'data/{cytokine}_{group}.tsv'
        tree = tree
        gene_key = config['gene_key']
    output:
        file_name = '{cytokine}_{group}_{geno}_hogwash_grouped'
        dir = "results/{cytokine}/results/"
    log:
        "log/{cytokine}_{group}_{geno}_hogwash_grouped.txt"
    resources:
        ncores = ncores
    scripts:
        "code/run_hogwash_grouped.R"

rule elastic_net:
    input: "combined_core.tsv"
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
