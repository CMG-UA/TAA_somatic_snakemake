import os

CONDITIONS = ["1", "2"]

# Define directories
REFDIR = os.getcwd()
#print(REFDIR)
sample_dir = REFDIR+"/data"

sample_names = []
sample_list = os.listdir(sample_dir)
for i in range(len(sample_list)):
    sample = sample_list[i]
    if sample.endswith("_1.fastq.gz"):
        samples = sample.split("_1.fastq")[0]
        sample_names.append(samples)
        #print(sample_names)

rule all:
    input:"results/01_multiqc/multiqc_report.html",
        "results/04_ATmultiqc/multiqc_report.html",
        expand("results/07_picard/marked_duplicates_{names}.bam", names=sample_names),
        expand("results/07_picard/marked_dup_metrics{names}.txt", names=sample_names)
        

rule fastqc: 
    input:
        "data/{names}_{con}.fastq.gz"
    output:
        result = directory("results/00_fastqc/{names}_{con}_fastqc/")
    log:
        "logs/fastqc_{names}_{con}.log"
    params:
        extra="-t 32"
    shell:
        """
        fastqc {params.extra} {input} --extract -o results/00_fastqc/ 2>> {log}
        """

rule multiqc:
    input:
        expand("results/00_fastqc/{names}_{con}_fastqc/", names=sample_names, con = CONDITIONS)
    output:
        "results/01_multiqc/multiqc_report.html",
        result = directory("results/01_multiqc/")  
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc results/00_fastqc/ -o {output.result} 2>> {log}
        """

rule Fastp:
    input:
        first = "data/{names}_1.fastq.gz",
        second = "data/{names}_2.fastq.gz"
    output:
        first = "results/02_fastp/{names}_1.fastq.gz",
        second = "results/02_fastp/{names}_2.fastq.gz",
        html = "results/02_fastp/{names}_fastp.html",
        json = "results/02_fastp/{names}_fastp.json"
    params:
        extra="-w 16"
    log:
        "logs/fastp_{names}.log"
    shell:
        """
        fastp {params.extra} -i {input.first} -I {input.second} -o {output.first} -O {output.second} -h {output.html} -j {output.json} --detect_adapter_for_pe 2>> {log}
        """
rule aftertrimming_fastqc: 
    input:
        "results/02_fastp/{names}_{con}.fastq.gz"
    output:
        result = directory("results/03_ATfastqc/{names}_{con}_fastqc/")
    log:
        "logs/ATfastqc_{names}_{con}.log"
    params:
        extra="-t 32"
    shell:
        """
        fastqc {params.extra} {input} --extract -o results/03_ATfastqc/ 2>> {log}
        """

rule aftertrimming_multiqc:
    input:
        expand("results/03_ATfastqc/{names}_{con}_fastqc/", names=sample_names, con = CONDITIONS)
    output:
        "results/04_ATmultiqc/multiqc_report.html",
        result = directory("results/04_ATmultiqc/")  
    log:
        "logs/ATmultiqc.log"
    shell:
        """
        multiqc results/03_ATfastqc/ -o {output.result} 2>> {log}
        """

rule BWA_indexing: 
    input: 
        "data/ref_data/hg38.fasta"
    output: 
        "data/ref_data/hg38.fasta.amb",
        "data/ref_data/hg38.fasta.ann",
        "data/ref_data/hg38.fasta.bwt",
        "data/ref_data/hg38.fasta.pac",
        "data/ref_data/hg38.fasta.sa"
    log: 
        "logs/bwa_index.log"
    shell: 
        """
        bwa index {input} 2>> {log}
        """

rule BWA_alignment:
    input: 
        "data/ref_data/hg38.fasta.amb",
        "data/ref_data/hg38.fasta.ann",
        "data/ref_data/hg38.fasta.bwt",
        "data/ref_data/hg38.fasta.pac",
        "data/ref_data/hg38.fasta.sa",
        first = "results/02_fastp/{names}_1.fastq.gz",
        second = "results/02_fastp/{names}_2.fastq.gz"
    output:
        "results/05_bwa/{names}.sam"
    log: 
        "logs/bwa_{names}.log"
    params:
        extra="-t 16"
    shell:
     """
     bwa mem {params.extra} data/ref_data/hg38.fasta {input.first} {input.second} > {output} 2>> {log}
     """

rule samtools: 
    input:
        "results/05_bwa/{names}.sam"
    output: 
        "results/06_samtools/{names}.bam"
    log:
        "logs/samtools_{names}.log"
    params:
        threads="-@ 32"
    shell:
        """
        samtools sort {params.threads} {input} -o {output} 2>> {log}
        samtools index {params.threads} {output} 2>> {log}
        """ 

rule Picard: 
    input:
        "results/06_samtools/{names}.bam"
    output:
        bam="results/07_picard/marked_duplicates_{names}.bam",
        txt="results/07_picard/marked_dup_metrics{names}.txt"
    log: 
        "logs/picard_{names}.log"
    params: 
        regex="--READ_NAME_REGEX null"
    shell: 
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.txt} {params.regex} 2>> {log}
        """