# Variant Calling Snakemake Pipeline

## Overview
This repository provides a Snakemake-based variant calling pipeline using GATK, BWA, Samtools, Trimmomatic, FastQC, and MultiQC. It automates the process from raw FASTQ files to high-confidence, filtered variants in VCF format.

## Workflow Steps
1. Read QC (FastQC)
2. Adapter trimming (Trimmomatic)
3. Alignment (BWA MEM)
4. BAM processing (sorting, indexing, read groups)
5. Variant calling (GATK HaplotypeCaller)
6. Variant filtering (GATK VariantFiltration)
7. MultiQC summary

## Quickstart

1. **Edit `config.yaml` with your sample FASTQ paths and reference genome.**

2. **Dry run to check pipeline:**
    ```bash
    sbatch dryrun.sh
    ```

3. **Run the full pipeline:**
    ```bash
    sbatch snakemake_run.sh
    ```
## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/en/latest/)
- SLURM

## Folder Structure

- `Snakefile` — Main workflow file
- `config.yaml` — Sample and reference paths
- `envs/` — Conda environment YAMLs for reproducibility
- `data/` — (Optional) Test/sample data folder (do not include large files!)
- `dryrun.sh` — Script for pipeline dry-run
- `snakemake_run.sh` — Script to run pipeline (edit for cluster if needed)

## Citation

If you use this pipeline, please cite GATK, Snakemake, and all component tools.
