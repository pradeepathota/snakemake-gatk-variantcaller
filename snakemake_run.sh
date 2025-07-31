#!/bin/bash
#SBATCH --job-name=snakemake_run
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -p general
#SBATCH -A <account id>

# Load conda
source <path_to_conda>/etc/profile.d/conda.sh
conda activate <your_env_name>

# Navigate to the workflow directory
cd <path_to_your_project_directory>

# Run Snakemake with real execution
snakemake --use-conda --conda-frontend conda --cores 4 --configfile config.yaml --latency-wait 60 --rerun-incomplete
