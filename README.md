# Snakemake GATK Variant Calling Pipeline

## Project Objective

To design and execute a scalable and reproducible variant calling pipeline using Snakemake and GATK that includes read trimming, alignment, sorting, indexing, and variant calling for human genomic data.

---
## Tools & Dependencies

- **Snakemake** (workflow management)
- **GATK** (variant calling)
- **BWA** (read alignment)
- **Trimmomatic** (read trimming)
- **Samtools** (file manipulation)
- **MultiQC** and **FastQC** (quality control)
- **Conda** (environment management)
- **HPC Cluster with SLURM** (job scheduling)

---

## 📁 Directory Structure

```bash
GATK/
├── Snakefile                # Snakemake workflow
├── config.yaml              # Configuration file with sample IDs and paths
├── envs/                    # Conda environment YAMLs for each tool
├── data/                    # Input FASTQ files
├── reference/               # Reference genome files (e.g., .fa, .dict, .fai) not included here
├── results/                 # Output directory with trimmed, aligned, and variant files
├── logs/                    # Logs from each Snakemake rule
├── snakemake_run.sh         # Submission script to run Snakemake on HPC
├── dryrun.sh                # SLURM script for dry-run
└── .gitignore               # To avoid pushing large or sensitive files
```

##  Workflow
The pipeline performs the following steps:

1. **Quality Control** with FastQC
2. **Trimming** adapters and low-quality reads using Trimmomatic
3. **Alignment** to the reference genome using BWA
4. **Sorting and Indexing** the BAM files using Samtools
5. **Variant Calling** using GATK HaplotypeCaller
6. **Generating MultiQC Reports** for consolidated QC

## pipeline execution on HPC
1. **Dry run**
A dry run in Snakemake is a way to test and visualize the workflow without actually executing any commands. It ensures everything is configured correctly before launching a potentially long and expensive run.
```
snakemake -n -p --use-conda --configfile config.yaml
```
2. **Actual run**
```
snakemake --use-conda --cores 4 --latency-wait 60 --rerun-incomplete --configfile config.yaml
```
## slurm submission
```
sbatch snakemake_run.sh
```
## Results

The pipeline produces the output files for each step but i have kept the key output from this pipeline located in the `results/` directory:

- `*.vcf`: Variant call format files containing identified SNPs and INDELs
- `multiqc_report.html`: Consolidated QC report across all samples

