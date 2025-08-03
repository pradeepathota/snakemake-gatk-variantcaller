#!/bin/bash
#SBATCH --job-name=snakemake_run
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH -A r00239
#SBATCH -p general
#SBATCH --mail-user=vathota@iu.edu
#SBATCH --mail-type=ALL 

# Load Conda and activate your base or Snakemake-specific environment
source /N/u/vathota/Quartz/miniconda3/etc/profile.d/conda.sh
conda activate bio-env

# Navigate to your project directory (update path if needed)
cd /N/u/vathota/Quartz/gatk

# Run Snakemake with real execution
snakemake --use-conda --conda-frontend conda --cores 8 --configfile config.yaml --latency-wait 60 --rerun-incomplete /N/scratch/vathota/gatk/results/SRR622461/variants.filtered.vcf
