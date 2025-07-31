#!/bin/bash
#SBATCH --job-name=snakemake_dryrun
#SBATCH --output=logs/dryrun_%j.out
#SBATCH --error=logs/dryrun_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -A <your_account_id>
#SBATCH -p general

# Load conda
source <path_to_conda>/etc/profile.d/conda.sh
conda activate <your_env_name>

# Navigate to the workflow directory
cd <path_to_your_project_directory>

# Perform Snakemake dry run
snakemake -n -p --use-conda --configfile config.yaml


