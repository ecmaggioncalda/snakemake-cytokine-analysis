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

#EVALUATE SIGNIFICANT ELIX ASSOCIATIONS ----
cytokine_cols <- grep("serum|stool",
                      colnames(pheno_cov),
                      value = TRUE)

variable_adjustment <- vector("list", length(cytokine_cols))
names(variable_adjustment) <- cytokine_cols

for(i in 1:length(cytokine_cols)){
  
  model <- paste0("summary(lm(",
                  cytokine_cols[i],
                  " ~ elix_weighted_update, data = pheno_cov))")
  
  out_model <- eval(str2lang(model))
  
  out_coefs <- out_model$coefficients
    
    if(out_model$r.squared > 0.10){

      variable_adjustment[[cytokine_cols[i]]] <- out_coefs[,1]

    }else{
    
      variable_adjustment[[cytokine_cols[i]]] <- NA
  }
}

variable_adjustment <- variable_adjustment[!is.na(variable_adjustment[cytokine_cols])]

print("Cytokines to be evaluated:")
paste0(cytokine_cols)
print("Cytokines to be additionally evaluated when adjusted by weighted elix score:")
paste0(names(variable_adjustment))

unadjusted_out <- c()

#i=2

for(i in 1:length(cytokine_cols)){
  
  unadjusted_out <- pheno_cov %>%
    select(genome_id,
           cytokine_cols[i]) %>%
    drop_na()
  
  write_tsv(unadjusted_out,
            file = paste0("../data/", cytokine_cols[i], "_raw.tsv"))
  
}

adjusted_out <- c()

#i=1

for(i in 1:length(variable_adjustment)){
  
  adjusted_out <- pheno_cov %>%
    mutate(adjusted_cytokine = eval(str2lang(names(variable_adjustment)[i])) - (elix_weighted_update*variable_adjustment[[i]][2] + variable_adjustment[[i]][1])) %>%
    select(genome_id,
           adjusted_cytokine) %>%
    drop_na()
  
  colnames(adjusted_out) <- c("genome_id",
                              names(variable_adjustment)[i])
  
  write_tsv(adjusted_out,
            file = paste0("../data/", names(variable_adjustment)[i], "_adjusted.tsv"))
  
}

dir.create(paste0("../results/", cytokine_cols))
dir.create(paste0("../results/", cytokine_cols, "/results"))
dir.create(paste0("../results/", cytokine_cols, "/plots"))