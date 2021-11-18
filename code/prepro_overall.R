#Script to build the individual phenotype files for each cytokine
#Determine which cytokines will have an additional branch using weighted elix score adjustment
#This is a snakemake compatible script based on the script 02_patient_covariates_analysis.R and patient_metadata_summary_markdown.Rmd
# 
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
#               cytokine = c("serum_TNFa",
#                            "serum_IL.6",
#                            "serum_IL.8"),
#               patient_metadata = c("data/patient_metadata.csv"),
#               group = c("raw",
#                         "adjusted"),
#               core = c("data/combined_core.tsv"),
#               pan = c("data/pan_geno.tsv"),
#               tree = c("data/cytokine_rooted_tree.tree"),
#               genome = c("core",
#                          "pan"),
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
#                            pheno = c(paste0("data/pheno/", config@cytokine[1][[1]], "/", config@group[1][[1]],".tsv"))),
#                  wildcards = c(cytokine = config@cytokine[1][[1]],
#                                genome = c("pan")),
#                  log = c(paste0("log/", config@cytokine[1][[1]],"_prepro_overall.txt")),
#                  resources = c(ncores = config@ncores[1][[1]]),
#                  params = c(core_path = config@geno[['core']],
#                             pan_path = config@geno[['pan']],
#                             tree = config@tree[1][[1]]))

source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno_cov <- readr::read_csv(file = snakemake@input[["patient_metadata"]])
#pheno_cov <- readr::read_csv(file = "./data/patient_metadata.csv")

#dir.create("data/pheno")

#CYTOKINES----
cytokine_cols <- snakemake@wildcards[["cytokine"]]
#cytokine_cols <- grep("serum|stool", colnames(pheno_cov), value = TRUE)
# cytokine_cols <- "testing"
# cytokine_cols <- "serum_IL.2Ra"

#cytokine_cols <- cytokine_cols[1]
paths <- paste0(c("data/pheno/",
                  "results/",
                  "benchmarks/",
                  "figures/",
                  "log/",
                  "data/mikropml/"),
                cytokine_cols)

for(i in 1:length(paths)){
  
  if(dir.exists(paths[i]) == FALSE){
    
    dir.create(paths[i])
    
  }else{
    
    print(paste0("directory ", paths[i], " already exists"))
    
  }
}

# new_dir1 <- paste0("data/pheno/", cytokine_cols)
# dir.create(new_dir1)
# 
# new_dir2 <- paste0("results/", cytokine_cols)
# dir.create(new_dir2)
# 
new_dir3 <- paste0("results/", cytokine_cols, "/runs")
dir.create(new_dir3)
# 
# new_dir4 <- paste0("benchmarks/", cytokine_cols)
# dir.create(new_dir4)
# 
# new_dir5 <- paste0("figures/", cytokine_cols)
# dir.create(new_dir5)
# 
# new_dir6 <- paste0("log/", cytokine_cols)
# dir.create(new_dir6)
# 
# new_dir7 <- paste0("data/mikropml/", cytokine_cols)
# dir.create(new_dir7)

#UNADJUSTED FILE ----
#for(i in 1:length(cytokine_cols)){
  
unadjusted_out <- pheno_cov %>%
  filter(clade != "cryptic") %>%
  select(genome_id,
         all_of(cytokine_cols),
         elix_weighted_update) %>%
  drop_na() %>%
  mutate(elix_weighted_update = NULL)
  
  write_tsv(unadjusted_out,
            file = paste0("data/pheno/", cytokine_cols, "/raw.tsv"))
  
#}

#EVALUATE SIGNIFICANT ELIX ASSOCIATIONS, GENERATE ADJUSTED ----
#for(i in 1:length(cytokine_cols)){
  
  model <- paste0("summary(lm(",
                  cytokine_cols,
                  " ~ elix_weighted_update, data = pheno_cov))")
  
  out_model <- eval(str2lang(model))
  
  out_coefs <- out_model$coefficients
  
  if(out_model$r.squared > 0.10){
    
    adjusted_out <- pheno_cov %>%
      filter(clade != "cryptic") %>%
      mutate(adjusted_cytokine = eval(str2lang(cytokine_cols)) - (elix_weighted_update*out_coefs[2,1] + out_coefs[1,1])) %>%
      #view()
      select(genome_id,
             adjusted_cytokine,
             elix_weighted_update) %>%
      drop_na() %>%
      mutate(elix_weighted_update = NULL) #%>%
      #view()
    
    colnames(adjusted_out) <- c("genome_id",
                                cytokine_cols)
    
    write_tsv(adjusted_out,
              file = paste0("data/pheno/", cytokine_cols, "/adjusted.tsv"))
  }
#}
