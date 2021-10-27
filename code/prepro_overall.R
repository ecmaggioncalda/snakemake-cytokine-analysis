#Script to build the individual phenotype files for each cytokine
#Determine which cytokines will have an additional branch using weighted elix score adjustment
#This is a snakemake compatible script based on the script 02_patient_covariates_analysis.R and patient_metadata_summary_markdown.Rmd

source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno_cov <- readr::read_csv(file = snakemake@input[["patient_metadata"]])
#pheno_cov <- readr::read_csv(file = "../../data/patient_metadata.csv")

#UNADJUSTED FILE ----
unadjusted_out <- pheno_cov %>%
  select(genome_id,
         snakemake@params[["cytokine"]]) %>%
  drop_na()

write_tsv(unadjusted_out,
          file = paste0("../data/", snakemake@params[["cytokine"]], "_raw.tsv"))

#EVALUATE SIGNIFICANT ELIX ASSOCIATIONS ----

model <- paste0("summary(lm(",
                snakemake@params[["cytokine"]],
                  " ~ elix_weighted_update, data = pheno_cov))")
  
out_model <- eval(str2lang(model))
  
out_coefs <- out_model$coefficients
    
  if(out_model$r.squared > 0.10){

    adjusted_out <- pheno_cov %>%
      mutate(adjusted_cytokine = eval(str2lang(snakemake@params[["cytokine"]])) - (elix_weighted_update*variable_adjustment[[i]][2] + variable_adjustment[[i]][1])) %>%
      select(genome_id,
             adjusted_cytokine) %>%
      drop_na()
    
    colnames(adjusted_out) <- c("genome_id",
                                snakemake@params[["cytokine"]])
    
    write_tsv(adjusted_out,
              file = paste0("../data/", snakemake@params[["cytokine"]], "_adjusted.tsv"))
  }

dir.create(paste0("../results/", snakemake@params[["cytokine"]]))
dir.create(paste0("../results/", snakemake@params[["cytokine"]], "/results"))
dir.create(paste0("../results/", snakemake@params[["cytokine"]], "/plots"))
