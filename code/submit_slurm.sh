#!/bin/bash

#SBATCH --job-name=cyokine_GWAS

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100MB
#SBATCH --time=96:00:00

#SBATCH --output=log/hpc/slurm-%j_%x.out

#SBATCH --account=esnitkin1
#SBATCH --partition=standard

#SBATCH --mail-user=emilycma@umich.edu
#SBATCH --mail-type=BEGIN,END

time snakemake --profile config/slurm --latency-wait 90 --configfile config/config.yml --show-failed-logs --verbose --keep-incomplete
