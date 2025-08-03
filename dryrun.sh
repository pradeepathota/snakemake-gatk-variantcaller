#!/bin/bash
#SBATCH --job-name=snakemake_dryrun
#SBATCH --output=logs/dryrun_%j.out
#SBATCH --error=logs/dryrun_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -A r00239
#SBATCH -p general

source /N/u/vathota/Quartz/miniconda3/etc/profile.d/conda.sh
conda activate bio-env

# Navigate to your project directory
cd /N/u/vathota/Quartz/GATK

# Perform Snakemake dry run
snakemake -n -p --use-conda --configfile config.yaml


