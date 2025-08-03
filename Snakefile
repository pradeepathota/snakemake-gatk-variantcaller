configfile: "config.yaml"

SAMPLES = config["samples"]
REF = config["ref"]

def outdir(path):
    return config["output_dir"] + "/" + path

rule all:
    input:
        expand(outdir("{sample}/variants.filtered.vcf"), sample=SAMPLES),
        outdir("multiqc_report.html")

rule fastqc:
    input:
        R1 = lambda wc: SAMPLES[wc.sample]["R1"],
        R2 = lambda wc: SAMPLES[wc.sample]["R2"]
    output:
        html1 = outdir("{sample}/{sample}_1_fastqc.html"),
        html2 = outdir("{sample}/{sample}_2_fastqc.html")
    conda: "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {config[output_dir]}/{wildcards.sample}
        fastqc {input.R1} {input.R2} -o {config[output_dir]}/{wildcards.sample}/
        """
rule trim_reads:
    input:
        R1 = lambda wc: SAMPLES[wc.sample]["R1"],
        R2 = lambda wc: SAMPLES[wc.sample]["R2"]
    output:
        R1_trimmed = outdir("{sample}/R1_trimmed.fastq.gz"),
        R2_trimmed = outdir("{sample}/R2_trimmed.fastq.gz")
    conda: "envs/trimmomatic.yaml"
    shell:
        """
        mkdir -p {config[output_dir]}/{wildcards.sample}
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
        R1 = outdir("{sample}/R1_trimmed.fastq.gz"),
        R2 = outdir("{sample}/R2_trimmed.fastq.gz")
    output:
        bam = outdir("{sample}/aligned.bam")
    conda: "envs/bwa.yaml"
    shell:
        """
        mkdir -p {config[output_dir]}/{wildcards.sample}
        bwa mem {REF} {input.R1} {input.R2} | samtools view -Sb - > {output.bam}
        """
rule sort_bam:
    input:
        bam = outdir("{sample}/aligned.bam")
    output:
        sorted_bam = outdir("{sample}/sorted.bam")
    conda: "envs/samtools.yaml"
    shell:
        "samtools sort {input.bam} -o {output.sorted_bam}"

rule index_bam:
    input:
        sorted_bam = outdir("{sample}/sorted.bam")
    output:
        bai = outdir("{sample}/sorted.bam.bai")
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input.sorted_bam}"
rule add_read_groups:
    input:
        bam=outdir("{sample}/sorted.bam")
    output:
        bam=outdir("{sample}/sorted_RG.bam")
    params:
        RGID="1",
        RGLB="lib1",
        RGPL="illumina",
        RGPU="unit1"
    conda:
        "envs/picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam} \
            RGID={params.RGID} \
            RGLB={params.RGLB} \
            RGPL={params.RGPL} \
            RGPU={params.RGPU} \
            RGSM={wildcards.sample}
        """
rule index_rg_bam:
    input:
        bam=outdir("{sample}/sorted_RG.bam")
    output:
        bai=outdir("{sample}/sorted_RG.bam.bai")
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"

rule call_variants:
    input:
        bam=outdir("{sample}/sorted_RG.bam"),
        bai=outdir("{sample}/sorted_RG.bam.bai"),
        ref=REF
    output:
        outdir("{sample}/variants.vcf")
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output}
        """
rule filter_variants:
    input:
        vcf=outdir("{sample}/variants.vcf"),
        ref=REF
    output:
        filtered_vcf=outdir("{sample}/variants.filtered.vcf")
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.filtered_vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum"
        """
rule multiqc:
    input:
        expand(outdir("{sample}/{sample}_1_fastqc.html"), sample=SAMPLES),
        expand(outdir("{sample}/{sample}_2_fastqc.html"), sample=SAMPLES)
    output:
        outdir("multiqc_report.html")
    conda: "envs/multiqc.yaml"
    shell:
        "multiqc {config[output_dir]} -o {config[output_dir]}"
