#Script to run treeWAS on each cytokine/genome combination
# setClass("config",
#          slots = c(ncores = "numeric",
#                    cytokine = "character",
#                    patient_metadata = "character",
#                    group = "character",
#                    core = "character",
#                    pan = "character",
#                    tree = "character",
#                    genome = "character",
#                    gene_key = "character",
#                    ml_methods = "character",
#                    kfold = "numeric",
#                    nseeds = "numeric"))
# 
# config <- new("config",
#               ncores = c(4),
#               nseeds = c(2),
#               kfold = c(5),
#               cytokine = c("serum_IL.2Ra"),
#               patient_metadata = c("data/patient_metadata.csv"),
#               group = c("raw",
#                         "adjusted"),
#               core = c("data/combined_core.tsv"),
#               pan = c("data/pan_geno.tsv"),
#               tree = c("data/cytokine_rooted_tree.tree"),
#               genome = c("pan",
#                          "core"),
#               gene_key = c("data/grouped_core_gene_key.tsv"))
# 
# setClass("snakemake", slots = c(input = "character",
#                                 output = "character",
#                                 wildcards = "character",
#                                 log = "character",
#                                 resources = "numeric",
#                                 params = "character"))
# snakemake <- new("snakemake",
#                  input = c(R = "code/prepro_overall.R",
#                            patient_metadata = config@patient_metadata[1][[1]],
#                            pheno = c(paste0("data/pheno/",
#                                             config@cytokine[1][[1]],
#                                             "/",
#                                             config@group[1][[1]],
#                                             ".tsv")),
#                            csv = "test_frame.csv",
#                            rds = paste0("data/mikropml/",
#                                         config@cytokine[1][[1]],
#                                         "/",
#                                         config@group[1][[1]],
#                                         ".",
#                                         config@genome[1][[1]],
#                                         ".dat_proc.Rds")),
#                  wildcards = c(cytokine = config@cytokine[1][[1]],
#                                genome = config@genome[1][[1]]),
#                  log = c(paste0("log/", config@cytokine[1][[1]],"_prepro_overall.txt")),
#                  resources = c(ncores = config@ncores[1][[1]]),
#                  params = c(core_path = config@core[1][[1]],
#                             pan_path = config@pan[1][[1]],
#                             tree = config@tree[1][[1]],
#                             file_name = "testing_name",
#                             dir = c(paste0("results/", config@cytokine[1][[1]]))))

source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(ape)
library(tidyverse)
library(treeWAS)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#USE MIKROPML RDS ----
rds <- readRDS(snakemake@input[["rds"]])$dat_transformed
# rds <- readRDS("data/mikropml/serum_IL.2Ra/raw.pan.dat_proc.Rds")$dat_transformed

#PHENO ----
paste0(snakemake@input[["pheno"]])
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")
#pheno <- readr::read_delim(file = "../data/serum_IL.6_raw.tsv",
#                           delim = "\t")

#GENO ----
geno <- rds[, !grepl(snakemake@wildcards[['cytokine']], colnames(rds))]

# if(snakemake@wildcards[['genome']] == "pan"){
#   
#   print(paste0("using pan genome path:", snakemake@params[['pan_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['pan_path']],
#                      row.names = 1)
#   
# }else{
#   
#   print(paste0("using core genome path:", snakemake@params[['core_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['core_path']],
#                      row.names = 1)
# }

# geno <- read.delim(file = "data/pan_geno.tsv",
#                    row.names = 1)
# geno <- read.delim(file = "../data/pan_matrix/pan_geno.tsv",
#                    row.names = 1)

#CLEAN GENO & PHENO ----
geno <- as.data.frame(geno)
rownames(geno) <- pheno$genome_id

# pheno <- pheno %>%
#   filter(genome_id %in% colnames(geno))

pheno <- deframe(pheno)

# TREE ----
tree <- read.tree(snakemake@params[["tree"]])
#tree <- read.tree("../../data/cytokine_rooted_tree.tree")

#PRINT STATEMENTS ----
print("dimension of geno, nrows (variants), ncol (genomes)")
dim(geno)

#geno_bad <- geno[,-1]

#if(sum(colnames(geno_bad) %in% names(pheno)) != ncol(geno_bad) | sum(colnames(geno_bad) %in% tree$tip.label) != ncol(geno_bad) | sum(names(pheno) %in% tree$tip.label) != ncol(geno_bad)){
#  stop("Mismatch in geno/pheno/tree contents")
#}

if(sum(rownames(geno) %in% names(pheno)) != nrow(geno) | sum(rownames(geno) %in% tree$tip.label) != nrow(geno) | sum(names(pheno) %in% tree$tip.label) != nrow(geno)){
  stop("Mismatch in geno/pheno/tree contents")
}

# TREEWAS ----
treeWAS.out <- treeWAS(snps = geno,
                       phen = pheno,
                       tree = tree,
                       filename.plot = snakemake@output[["plot"]])

save(treeWAS.out, file = snakemake@output[["rdata"]])