#Script to run treeWAS on each cytokine/genome combination
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(ape)
library(tidyverse)
library(treeWAS)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#PHENO ----
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")
#pheno <- readr::read_delim(file = "../data/serum_IL.6_raw.tsv",
#                           delim = "\t")

#GENO ----
geno <- read.delim(file = snakemake@params[["geno"]],
                   row.names = 1)
#geno <- read.delim(file = "../../data/pan_matrix/pan_geno.tsv",
#                   row.names = 1)

#CLEAN GENO & PHENO ----
pheno <- pheno %>%
  filter(genome_id %in% colnames(geno))

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

if(sum(colnames(geno) %in% names(pheno)) != ncol(geno) | sum(colnames(geno) %in% tree$tip.label) != ncol(geno) | sum(names(pheno) %in% tree$tip.label) != ncol(geno)){
  stop("Mismatch in geno/pheno/tree contents")
}

# TREEWAS ----
treeWAS.out <- treeWAS(snps = t(geno),
                       phen = pheno,
                       tree = tree,
                       filename.plot = snakemake@output[["plot"]])

save(treeWAS.out, file = snakemake@output[["rdata"]])