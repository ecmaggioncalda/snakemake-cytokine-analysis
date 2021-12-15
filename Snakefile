configfile: 'config/config.yml'

cytokine = config['cytokine']
tree = config['tree']
genome = config['genome']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule all:
    input:
        expand("aggregated/{cytokine}.runs.csv", cytokine=cytokine)

checkpoint prepro_overall:
    input:
        R = "code/prepro_overall.R",
        patient_metadata = config['patient_metadata']
    output:
        dat_dir = directory("data/pheno/{cytokine}")
    log:
        "log/{cytokine}/prepro_overall.txt"
    resources:
        ncores = ncores
    script:
        "code/prepro_overall.R"

rule generate_mikropml_df:
    input:
        R = "code/generate_mikropml_df.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv"
    output:
        file_name = "data/mikropml/{cytokine}/{group}.{genome}.csv"
    params:
        snp_path = config['snp'],
        pan_path = config['pan'],
        gene_path = config['gene']
    log:
        "log/{cytokine}/{group}.{genome}.generate_mikropml_df.txt"
    resources:
        ncores = ncores
    script:
        "code/generate_mikropml_df.R"

include: "mikropml.smk"

rule run_treeWAS:
    input:
        R = "code/run_treeWAS.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv",
        rds = rules.preprocess_data.output.rds
    output:
        rdata = 'results/{cytokine}/{group}.{genome}.treeWAS.RData',
        plot = 'results/{cytokine}/{group}.{genome}.treeWAS.pdf'
    params:
        tree = tree
    log:
        "log/{cytokine}/{group}.{genome}.treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

rule run_hogwash_ungrouped:
    input:
        R = "code/run_hogwash_ungrouped.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv",
        rds = rules.preprocess_data.output.rds
    output:
        rdata = "results/{cytokine}/hogwash_continuous_{group}.{genome}.ungrouped.rda",
        plot = "results/{cytokine}/hogwash_continuous_{group}.{genome}.ungrouped.pdf"
    params:
        tree = tree,
        file_name = '{group}.{genome}.ungrouped',
        dir = "results/{cytokine}"
    log:
        "log/{cytokine}/{group}.{genome}.hogwash.ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

# rule run_hogwash_grouped:
#     input:
#         R = "code/run_hogwash_grouped.R",
#         pheno = "data/pheno/{cytokine}/{group}.tsv",
#         rds = rules.preprocess_data.output.rds
#     output:
#         rdata = protected("results/{cytokine}/hogwash_continuous_{group}.{genome}.grouped.rda"),
#         plot = protected("results/{cytokine}/hogwash_continuous_{group}.{genome}.grouped.pdf")
#     params:
#         tree = tree,
#         file_name = '{group}.{genome}.grouped',
#         dir = "results/{cytokine}",
#         gene_key = config['gene_key']
#     wildcard_constraints:
#         genome = "snp"
#     log:
#         "log/{cytokine}/{group}.{genome}.hogwash.grouped.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_hogwash_grouped.R"

def aggregate_input1(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.{genome}.treeWAS.RData',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

def aggregate_input2(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/hogwash_continuous_{group}.{genome}.ungrouped.rda',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

# def aggregate_input3(wildcards):
#     checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
#     return expand('results/{cytokine}/hogwash_continuous_{group}.{genome}.grouped.rda',
#         cytokine=wildcards.cytokine,
#         group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
#         genome = "snp")

def aggregate_input4(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.{genome}.report.md',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

# if 'snp' in genome:
#     finish_list = [aggregate_input1, aggregate_input2, aggregate_input3, aggregate_input4]
# else:
#     finish_list = [aggregate_input1, aggregate_input2, aggregate_input4]

rule finish_test:
    input:
        # finish_list
        aggregate_input1,
        aggregate_input2,
        aggregate_input4
    output:
        "aggregated/{cytokine}.runs.csv"
    log:
        "log/{cytokine}/finish.txt"
    script:
        "code/assemble_files.py"
