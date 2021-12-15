#Script to run hogwash on each cytokine/genome combination
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(ape)
library(tidyverse)
library(hogwash)
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

#CLEAN GENO & PHENO ----
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$genome_id

geno <- as.data.frame(geno)
rownames(geno) <- pheno$genome_id

pheno <- pheno[, -1, drop = FALSE]

# TREE ----
tree <- read.tree(snakemake@params[["tree"]])
# tree <- read.tree("../../data/cytokine_rooted_tree.tree")

#ORDER DATA ----
reorder_pheno <- match(tree$tip.label, rownames(pheno))

pheno <- as.matrix(pheno[reorder_pheno, , drop = FALSE])

reorder_geno <- match(tree$tip.label, rownames(geno))

geno <- as.matrix(geno[reorder_geno,  , drop = FALSE])

if(all(rownames(pheno) == tree$tip.label) | all(rownames(geno) == tree$tip.label)){
  
  print("Geno/Pheno/Tree contents are in the correct order")
  
}else{
  
  stop("Mismatch in geno/pheno/tree contents")
  
}

'Number of samples and variants'
dim(geno)

#RUN HOGWASH ----
hogwash(pheno = pheno, 
        geno = geno, 
        tree = tree, 
        file_name = snakemake@params[["file_name"]], 
        dir = snakemake@params[["dir"]])
