configfile: "config.yaml"

SAMPLES = config["samples"]
REF = config["ref"]

rule all:
    input:
        expand("results/{sample}/variants.vcf", sample=SAMPLES),
        "results/multiqc_report.html"

rule fastqc:
    input:
        R1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        R2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        "results/{sample}/{sample}_1_fastqc.html",
        "results/{sample}/{sample}_2_fastqc.html"
    conda: "envs/fastqc.yaml"
    shell:
        "fastqc {input.R1} {input.R2} -o results/{wildcards.sample}/"


rule trim_reads:
    input:
        R1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        R2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        R1_trimmed = "results/{sample}/R1_trimmed.fastq.gz",
        R2_trimmed = "results/{sample}/R2_trimmed.fastq.gz"
    conda: "envs/trimmomatic.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}
        trimmomatic PE {input.R1} {input.R2} \
        {output.R1_trimmed} /dev/null \
        {output.R2_trimmed} /dev/null \
        SLIDINGWINDOW:4:20 MINLEN:36
        """

rule bwa_index:
    input: REF
    output: REF + ".bwt"
    conda: "envs/bwa.yaml"
    shell:
        "bwa index {input}"

rule align:
    input:
        idx = REF + ".bwt",
        R1 = "results/{sample}/R1_trimmed.fastq.gz",
        R2 = "results/{sample}/R2_trimmed.fastq.gz"
    output: "results/{sample}/aligned.bam"
    conda: "envs/bwa.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}
        bwa mem {REF} {input.R1} {input.R2} | samtools view -Sb - > {output}
        """

rule sort_bam:
    input: "results/{sample}/aligned.bam"
    output: "results/{sample}/sorted.bam"
    conda: "envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule index_bam:
    input: "results/{sample}/sorted.bam"
    output: "results/{sample}/sorted.bam.bai"
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule call_variants:
    input:
        bam = "results/{sample}/sorted.bam",
        bai = "results/{sample}/sorted.bam.bai",
        ref = REF
    output: "results/{sample}/variants.vcf"
    conda: "envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -O {output}
        """

rule multiqc:
    input:
        expand("results/{sample}/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/{sample}/{sample}_2_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report.html"
    conda: "envs/multiqc.yaml"
    shell:
        "multiqc results/ -o results/"
