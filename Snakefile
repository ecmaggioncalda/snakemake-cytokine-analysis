configfile: 'config/config.yml'

cytokine = config['cytokine']
tree = config['tree']
#group = ["raw", "adjusted"]
genome = ["pan", "core"]
model = ["hogwash.ungrouped.rda", "treeWAS.RData"]

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

# wildcard_constraints:
#     #group="\w",
#     genome="\w"
#
# ruleorder: prepro_overall > generate_mikropml_df

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
        "log/{cytokine}_prepro_overall.txt"
    resources:
        ncores = ncores
    script:
        "code/prepro_overall.R"

rule run_treeWAS:
    input:
        R = "code/run_treeWAS.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv"
    output:
        rdata = 'results/{cytokine}/{group}.{genome}.treeWAS.RData',
        plot = 'results/{cytokine}/{group}.{genome}.treeWAS.pdf'
    params:
        core_path = config['core'],
        pan_path = config['pan'],
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
        pheno = "data/pheno/{cytokine}/{group}.tsv"
    output:
        rdata = "results/{cytokine}/{group}.{genome}.hogwash.ungrouped.rda",
        plot = "results/{cytokine}/{group}.{genome}.hogwash.ungrouped.pdf"
    params:
        core_path = config['core'],
        pan_path = config['pan'],
        tree = tree,
        file_name = '{group}.{genome}.hogwash.ungrouped',
        dir = "results/{cytokine}"
    log:
        "log/{cytokine}/{group}.{genome}.hogwash.ungrouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_ungrouped.R"

rule run_hogwash_grouped:
    input:
        R = "code/run_hogwash_grouped.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv"
    output:
        rdata = "results/{cytokine}/{group}.core.hogwash.grouped.rda",
        plot = "results/{cytokine}/{group}.core.hogwash.grouped.pdf"
    params:
        core_path = config['core'],
        tree = tree,
        file_name = '{group}.core.hogwash.grouped',
        dir = "results/{cytokine}",
        gene_key = config['gene_key']
    log:
        "log/{cytokine}/{group}.core.hogwash.grouped.txt"
    resources:
        ncores = ncores
    script:
        "code/run_hogwash_grouped.R"

rule generate_mikropml_df:
    input:
        R = "code/generate_mikropml_df.R",
        pheno = "data/pheno/{cytokine}/{group}.tsv"
    output:
        file_name = "data/mikropml/{cytokine}/{group}.{genome}.csv"
    params:
        core_path = config['core'],
        pan_path = config['pan'],
    log:
        "log/{cytokine}/{group}.{genome}.generate_mikropml_df.txt"
    resources:
        ncores = ncores
    script:
        "code/generate_mikropml_df.R"

include: "rules/mikropml.smk"

def aggregate_input1(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.{genome}.{model}',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome,
        model=model)

def aggregate_input2(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.core.hogwash.grouped.rda',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group)

# def aggregate_input3(wildcards):
#     checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
#     return expand('data/mikropml/{cytokine}/{group}.{genome}.csv',
#         cytokine=wildcards.cytokine,
#         group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
#         genome=genome)

def aggregate_input4(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.{genome}.report.md',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

rule finish_test:
    input:
        aggregate_input1,
        aggregate_input2,
        # aggregate_input3,
        aggregate_input4
    output:
        "aggregated/{cytokine}.runs.csv"
    log:
        "log/{cytokine}_finish.txt"
    shell:
        "cat {input} > {output}"
    # "   import os       "
    # "   os.listdir('./results/{wildcards.cytokine}') > file_list        "
    # "   import csv      "
    # "   with open('aggregated/{cytokine}.runs.csv', 'wt', newline= '') as outfile:      "
    # "       writer = csv.writer(outfile, delimiter = '\t')      "
    # "       writer.writerows()     "
    #system
         # I should modify this so that the shell output is something that is actually useful, right now it is non-human readable but functions for the workflow
         # '''
         # import csv
         #     with open('aggregated/{cytokine}.tsv', 'wt') as out_file:
         #         tsv_writer = csv.writer(out_file, delimiter='\t')
         #         tsv_writer.writerows({input})
         # '''
#

#include: "rules/summarize.smk", this should probably actually be a directory up instead of below, but I didn't think it through this way initially
#alternatively could create a sub-workflow snakefile that just depends on the outputs, need to figure out which is preferable workflow-wise
#
#rule summarize_cytokine_results:
#    input:
#    output:
#    log:
#    params:
#    resources:
#    script:
#
#probably also want to mark the outputs of hogwash as protected since they take so long to run?
#
