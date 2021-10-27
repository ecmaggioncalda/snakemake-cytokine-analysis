#Script to run hogwash on each cytokine/genome combination
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(ape)
library(tidyverse)
library(hogwash)
library(doFuture)

#PHENO ----
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")
#pheno <- read.delim(file = "../data/serum_IL.6_raw.tsv",
#                    row.names = 1)

#GENO ----
geno <- read.delim(file = snakemake@params[["geno"]],
                   row.names = 1)
#geno <- read.delim(file = "../../data/pan_matrix/pan_geno.tsv",
#                   row.names = 1)

#CLEAN GENO & PHENO ----
pheno <- pheno[rownames(pheno) %in% colnames(geno), , drop = FALSE]

# TREE ----
tree <- read.tree(snakemake@params[["tree"]])
#tree <- read.tree("../../data/cytokine_rooted_tree.tree")

#ORDER DATA ----
reorder_pheno <- match(tree$tip.label, rownames(pheno))

pheno <- pheno[reorder_pheno, , drop = FALSE]

reorder_geno <- match(tree$tip.label, rownames(t(geno)))

geno <- t(geno[, reorder_geno , drop = FALSE])

if(all(rownames(pheno) == tree$tip.label) | all(rownames(geno) == tree$tip.label) == FALSE){
  stop("Mismatch in geno/pheno/tree contents")
}

'Number of samples and variants'
dim(geno)

#LOAD GENE KEY ----
gene_key <- snakemake@params[["gene_key"]]

#RUN HOGWASH ----
hogwash(pheno = pheno, 
        geno = geno, 
        tree = tree, 
        file_name = snakemake@output[["file_name"]],
        dir = snakemake@params[["dir"]],
        group_genotype_key = gene_key,
        grouping_method = "post-ar",
        test='continuous')