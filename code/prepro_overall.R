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
pheno_cov <- pheno_cov %>%
  drop_na(elix_weighted_update) %>%
  filter(clade != "cryptic")

# pheno_cov <- readr::read_csv(file = "./data/patient_metadata.csv")

#CYTOKINES----
cytokine_cols <- snakemake@wildcards[["cytokine"]]
# cytokine_cols <- grep("serum|stool", colnames(pheno_cov), value = TRUE)
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

new_dir3 <- paste0("results/", cytokine_cols, "/runs")

if(dir.exists(new_dir3) == FALSE){
  
  dir.create(new_dir3)
    
  }else{
    
    print(paste0("directory ", new_dir3, " already exists"))
    
  }



#UNADJUSTED FILE ----
unadjusted_out <- pheno_cov %>%
  filter(clade != "cryptic") %>%
  select(genome_id,
         all_of(cytokine_cols),
         elix_weighted_update) %>%
  drop_na() %>%
  mutate(elix_weighted_update = NULL)
  
  write_tsv(unadjusted_out,
            file = paste0("data/pheno/", cytokine_cols, "/raw.tsv"))

#EVALUATE SIGNIFICANT ELIX ASSOCIATIONS, GENERATE ADJUSTED ----
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
