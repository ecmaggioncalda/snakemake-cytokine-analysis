configfile: 'config/config.yml'

wildcard_constraints:
#    cytokine='\w+',
    group='\w+',
    genome='\w+',
    model='\w+'

cytokine = config['cytokine']

tree = config['tree']
core = config['core']
pan = config['pan']
ncores = config['ncores']

#testing vars
genome = ["pan"]
#genome = ["pan", "core"]

group = ["raw"]
#group = ["raw", "adjusted"]

model = ["treeWAS"]
#model = ["hogwash_ungrouped", "treeWAS.RData"]

rule all:
    input:
        expand("aggregated/{cytokine}.tsv", cytokine=cytokine)
        #"results/all_runs.tsv"
        #expand('results/{cytokine_file}.pan.treeWAS.RData')
        #"aggregated/{cytokine}.txt"
        #expand("results/{cytokine}/{group}.{genome}.{model}.RData", cytokine=cytokine, group=group, genome=genome, model=model)
        #expand("results/{cytokine}/results/{group}_core_hogwash_grouped", cytokine = cytokine, group = group, allow_missing=True),
        #expand("mikropml-snakemake-workflow/data/{cytokine}_{group}_{genome}.csv", cytokine = cytokine, group = group, genome = genome, allow_missing=True),
        #expand("mikropml-snakemake-workflow/{cytokine}_{group}_{genome}_report.md", cytokine = cytokine, group = group, genome = genome, allow_missing=True)

checkpoint prepro_overall:
#rule prepro_overall:
    input:
        R = "code/prepro_overall.R",
        patient_metadata = config['patient_metadata']
    output:
        dat_dir = directory("data/{cytokine}"),
        #res_dir = directory("results/{cytokine}")
        #dir = directory("results/{cytokine}")
        #dir = directory("data/pheno")
        #files = "data/pheno/{cytokine_file}.tsv"
    #params:
        #cytokine = config['cytokine']
    #log:
    #    "log/prepro_overall.txt"
    resources:
        ncores = ncores
    script:
        "code/prepro_overall.R"

rule run_treeWAS_pan:
    input:
        R = "code/run_treeWAS.R",
        pheno = "data/{cytokine}/{group}.tsv"
        #pheno = unpack(aggregate_input)
        #pheno = aggregate_input
        #pheno = cytokine_files
    output:
        #multiext('results/{cytokine}/{group}.pan.treeWAS', ".RData", ".pdf")
        rdata = 'results/{cytokine}/{group}.pan.treeWAS.RData',
        plot = 'results/{cytokine}/{group}.pan.treeWAS.pdf'
    params:
        geno = pan,
        tree = tree
    log:
        "log/{cytokine}/{group}.pan.treeWAS.txt"
    resources:
        ncores = ncores
    script:
        "code/run_treeWAS.R"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{cytokine}/{group}.tsv',
        cytokine=wildcards.cytokine,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group)
    #return cytokine_file = glob_wildcards(os.path.join(checkpoint_path, "{cytokine_file}.tsv"))
    #cytokine_file = glob_wildcards("data/pheno/{cytokine_file}.tsv")
    #return expand("data/pheno/{cytokine_file}.tsv", cytokine_file=cytokine_file)

rule finish_test:
    input:
        aggregate_input
        #expand('results/{cytokine_file}.pan.treeWAS.RData')
        #expand(rules.run_treeWAS_pan.output.rdata, cytokine_file=cytokine_file)
        #expand('results/{cytokine_file}.pan.treeWAS.RData', cytokine_file=glob_wildcards("data/pheno/{cytokine_file}.tsv"))
    output:
        "aggregated/{cytokine}.tsv"
    shell:
        "cat {input} > {output}"
         # '''
         # import csv
         #     with open('aggregated/{cytokine}.tsv', 'wt') as out_file:
         #         tsv_writer = csv.writer(out_file, delimiter='\t')
         #         tsv_writer.writerows({input})
         # '''

# def aggregate_input(wildcards):
#     cytokine_files = glob_wildcards("data/pheno/{cytokine_file}.tsv")
#     return expand("data/pheno/{cytokine_file}.tsv", cytokine_file=cytokine_files)
#     #return expand("results/{cytokine_file}.pan.treeWAS.RData",
#     #cytokine_file=glob_wildcards(checkpoints.prepro_overall.get(**wildcards).output[0]))
#     #checkpoint_output = glob_wildcards(checkpoints.prepro_overall.get(**wildcards).output[0])
#     #return expand("results/{cytokine_file}.pan.treeWAS.RData",
#     #cytokine_file=checkpoint_output)
#      #cytokine_file=glob_wildcards(os.path.join(checkpoint_output, "{cytokine_file}.tsv")).cytokine_file)

# rule aggregate:
#     input:
#         aggregate_input
#     output:
#         combined = "aggregated/cytokine_files.txt"
#     shell:
#         "cat {input} > {output.combined}"
#         #"import csv
#         #with open('aggregated/cytokine_files.tsv', 'w') as out_file:
#         #    tsv_writer = csv.writer(out_file, delimiter='\t')
#         #    tsv_writer.writerows({input})"

# rule run_hogwash_ungrouped_pan:
#     input:
#         R = "code/run_hogwash_ungrouped.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         "results/{cytokine}/results/{group}_pan_hogwash_ungrouped"
#     params:
#         file_name = '{cytokine}_pan_hogwash_ungrouped',
#         geno = pan,
#         tree = tree,
#         dir = "results/{cytokine}/results/"
#     log:
#         "log/{cytokine}/{group}_pan_hogwash_ungrouped.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_hogwash_ungrouped.R"
#
# # rule run_treeWAS_core:
#     input:
#         R = "code/run_treeWAS.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         rdata = 'results/{cytokine}/results/{group}_core_treeWAS.RData',
#         plot = 'results/{cytokine}/plots/{group}_core_treeWAS.pdf'
#     params:
#         geno = core,
#         tree = tree
#     log:
#         "log/{cytokine}/{group}_core_treeWAS.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_treeWAS.R"
#
# rule run_hogwash_ungrouped_core:
#     input:
#         R = "code/run_hogwash_ungrouped.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         "results/{cytokine}/results/{group}_core_hogwash_ungrouped"
#     params:
#         file_name = '{group}_core_hogwash_ungrouped',
#         geno = core,
#         tree = tree,
#         dir = "results/{cytokine}/results/"
#     log:
#         "log/{cytokine}/{group}_core_hogwash_ungrouped.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_hogwash_ungrouped.R"
#
# rule run_hogwash_grouped:
#     input:
#         R = "code/run_hogwash_grouped.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         "results/{cytokine}/results/{group}_core_hogwash_grouped"
#     params:
#         file_name = '{cytokine}_core_hogwash_grouped',
#         geno = core,
#         tree = tree,
#         dir = "results/{cytokine}/results/",
#         gene_key = config['gene_key']
#     log:
#         "log/{cytokine}/{group}_core_hogwash_grouped.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/run_hogwash_grouped.R"
#
# rule generate_mikropml_df_core:
#     input:
#         R = "code/generate_mikropml_df.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         "mikropml-snakemake-workflow/data/{cytokine}_{group}_core.csv"
#     params:
#         geno = core
#     log:
#         "log/{cytokine}/{group}_generate_mikropml_df_core.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/generate_mikropml_df.R"
#
# rule generate_mikropml_df_pan:
#     input:
#         R = "code/generate_mikropml_df.R",
#         pheno = 'data/{cytokine}_{group}.tsv'
#     output:
#         "mikropml-snakemake-workflow/data/{cytokine}_{group}_pan.csv"
#     params:
#         geno = pan
#     log:
#         "log/{cytokine}/{group}_generate_mikropml_df_pan.txt"
#     resources:
#         ncores = ncores
#     script:
#         "code/generate_mikropml_df.R"
#
#include: "mikropml-snakemake-workflow/Snakefile"
#
#rule summarize_cytokine_results:
#    input:
#    output:
#    log:
#    params:
#    resources:
#    script:
#
###Come back to this, need a checkpoint to trim the initial DAG since not all cytokines have an adjusted file
###Documentation doesn't make it clear whether it is necessary for all wildcard combinations of the output files to be present for the workflow to continue, so it might be ok, the checkpoint information seems more specific to the generation of files as opposed to the absence
#checkpoint adjustments:
#    input:
#        R = "code/verify_adjustment_files.R"
#        "data/{cytokine}.tsv"
#    output:
#
#    script:
#        "code/verify_adjustment_files.R"
#
#probably also want to mark the outputs of hogwash as protected since they take so long to run?
