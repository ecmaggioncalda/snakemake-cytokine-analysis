#summary file
#Script to run hogwash on each cytokine/genome combination
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(ape)
library(tidyverse)
library(doFuture)

#Read in results ----
filenames <- as.list(list.files(path= snakemake@param[["dir"]],
                                pattern=".rda"))

for (i in 1:length(filenames)){
  
  load(file = paste(snakemake@param[["dir"]],
                  filenames[[i]],
                  sep=""))
  
  a <- load(file = paste(snakemake@param[["dir"]],
                       filenames[[i]],
                       sep=""))
  
  assign(paste(filenames[[i]],get(a)))
}

