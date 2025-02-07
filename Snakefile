import os

CONDITIONS = ["1", "2"]
TYPE = ["tumor", "blood"]

# Define directories
REFDIR = os.getcwd()
#print(REFDIR)
sample_dir = REFDIR+"/data"

sample_names = []
sample_list = os.listdir(sample_dir)
for i in range(len(sample_list)):
    sample = sample_list[i]
    if sample.endswith("_tumor_1.fastq.gz"):
        samples = sample.split("_tumor_1.fastq")[0]
        sample_names.append(samples)
        #print(sample_names)

rule all:
    input:"results/01_multiqc/multiqc_report.html",
        "results/04_ATmultiqc/multiqc_report.html",
        expand("results/09_variant_calling/mutect/{names}_somatic.vcf.gz", names=sample_names),
        expand("results/09_variant_calling/varscan/{names}.vcf", names=sample_names),
        expand("results/09_variant_calling/somaticsniper/{names}.vcf", names=sample_names),
        expand("results/09_variant_calling/MuSE/{names}.vcf", names=sample_names)
        

rule fastqc: 
    input:
        "data/{names}_{type}_{con}.fastq.gz"
    output:
        result = directory("results/00_fastqc/{names}_{type}_{con}_fastqc/")
    log:
        "logs/fastqc_{names}_{type}_{con}.log"
    params:
        extra="-t 64"
    shell:
        """
        fastqc {params.extra} {input} --extract -o results/00_fastqc/ 2>> {log}
        """

rule multiqc:
    input:
        expand("results/00_fastqc/{names}_{type}_{con}_fastqc/", names=sample_names, type=TYPE, con = CONDITIONS)
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
        first = "data/{names}_{type}_1.fastq.gz",
        second = "data/{names}_{type}_2.fastq.gz"
    output:
        first = "results/02_fastp/{names}_{type}_1.fastq.gz",
        second = "results/02_fastp/{names}_{type}_2.fastq.gz",
        html = "results/02_fastp/{names}_{type}_fastp.html",
        json = "results/02_fastp/{names}_{type}_fastp.json"
    params:
        extra="-w 16"
    log:
        "logs/fastp_{names}_{type}.log"
    shell:
        """
        fastp {params.extra} -i {input.first} -I {input.second} -o {output.first} -O {output.second} -h {output.html} -j {output.json} --detect_adapter_for_pe 2>> {log}
        """
rule aftertrimming_fastqc: 
    input:
        "results/02_fastp/{names}_{type}_{con}.fastq.gz"
    output:
        result = directory("results/03_ATfastqc/{names}_{type}_{con}_fastqc/")
    log:
        "logs/ATfastqc_{names}_{type}_{con}.log"
    params:
        extra="-t 64"
    shell:
        """
        fastqc {params.extra} {input} --extract -o results/03_ATfastqc/ 2>> {log}
        """

rule aftertrimming_multiqc:
    input:
        expand("results/03_ATfastqc/{names}_{type}_{con}_fastqc/", names=sample_names, type=TYPE, con = CONDITIONS)
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
        first = "results/02_fastp/{names}_{type}_1.fastq.gz",
        second = "results/02_fastp/{names}_{type}_2.fastq.gz"
    output:
        "results/05_bwa/{names}_{type}.sam"
    log: 
        "logs/bwa_{names}_{type}.log"
    params:
        extra="-t 64"
    shell:
     """
     bwa mem {params.extra} data/ref_data/hg38.fasta {input.first} {input.second} > {output} 2>> {log}
     """

rule samtools: 
    input:
        "results/05_bwa/{names}_{type}.sam"
    output: 
        "results/06_samtools/{names}_{type}.bam"
    log:
        "logs/samtools_{names}_{type}.log"
    params:
        threads="-@ 64"
    shell:
        """
        samtools sort {params.threads} {input} -o {output} 2>> {log}
        samtools index {params.threads} {output} 2>> {log}
        """ 
rule AddReadGroups:
    input:
        "results/06_samtools/{names}_{type}.bam"
    output:
        first="results/06_samtools/{names}_{type}_rg.bam",
        second="results/06_samtools/{names}_{type}_rg_sorted.bam"
    log:
        "logs/add_read_groups_{names}_{type}.log"
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input} \
            O={output.first} \
            RGID={wildcards.names}_{wildcards.type} \
            RGLB=lib1 \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM={wildcards.names}_{wildcards.type} \
            2>> {log}

            samtools sort -o {output.second} {output.first} 2>> {log}
        """

rule Picard: 
    input:
        "results/06_samtools/{names}_{type}_rg_sorted.bam"
    output:
        bam="results/07_picard/marked_duplicates_{names}_{type}.bam",
        txt="results/07_picard/marked_dup_metrics_{names}_{type}.txt"
    log: 
        "logs/picard_{names}_{type}.log"
    params:
        regex="--READ_NAME_REGEX null",
        resources="--MAX_RECORDS_IN_RAM 2000000 --COMPRESSION_LEVEL 5"
    shell: 
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.txt} {params.regex} {params.resources} 2>> {log}
        """
rule mpileup:
    input: 
        "results/07_picard/marked_duplicates_{names}_{type}.bam"
    output:
        "results/08_mpileup/marked_duplicates_{names}_{type}.pileup"
    log: 
        "logs/mpileup_{names}_{type}.log"
    params: 
        ref="-f data/ref_data/hg38.fasta"
    shell: 
        """
        samtools mpileup {params.ref} {input} > {output} 2>> {log}
        """

rule Varscan: 
    input: 
        disease="results/08_mpileup/marked_duplicates_{names}_tumor.pileup",
        blood="results/08_mpileup/marked_duplicates_{names}_blood.pileup"
    output: 
        snp="results/09_variant_calling/varscan/{names}.snp",
        indel="results/09_variant_calling/varscan/{names}.indel"
    log: 
        "logs/varscan_{names}.log"
    #params: 
    shell: 
        """
        varscan somatic {input.blood} {input.disease} --output-snp {output.snp} --output-indel {output.indel}
        """

rule snp_to_vcf:
    input:
        snp="results/09_variant_calling/varscan/{names}.snp"
    output:
        "results/09_variant_calling/varscan/{names}.vcf"
    log:
        "logs/convert_snp_vcf_{names}.log"
    shell:
        """
        python3 scripts/vs_format_converter.py {input.snp} >> {output} 2>> {log}
        """

rule fasta_dict: 
    input:
        ref="data/ref_data/hg38.fasta"
    output:
        "data/ref_data/hg38.dict"
    log:
        "logs/fasta_dict.log"
    shell:
        """
        gatk CreateSequenceDictionary -R {input.ref} 2>> {log}
        """
        
rule index:
    input:
        "results/07_picard/marked_duplicates_{names}_{type}.bam"
    output: 
        "results/07_picard/marked_duplicates_{names}_{type}.bam.bai"
    log: 
        "logs/index_dup_{names}_{type}.log"
    params: 
        "-@ 32"
    shell: 
        "samtools index {params} {input}"

rule Mutect: 
    input: 
        disease="results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="results/07_picard/marked_duplicates_{names}_blood.bam",
        index_d="results/07_picard/marked_duplicates_{names}_tumor.bam.bai",
        index_b="results/07_picard/marked_duplicates_{names}_blood.bam.bai",
        ref_dict="data/ref_data/hg38.dict"
    output: 
        somatic="results/09_variant_calling/mutect/{names}_somatic.vcf.gz"
    log: 
        "logs/mutect_{names}.log"
    params: 
        ref="-R data/ref_data/hg38.fasta"
    shell: 
        """
        gatk Mutect2 \
        {params.ref} \
        -I {input.disease}\
        -I {input.blood} \
        -normal {wildcards.names}_blood \
        -O {output.somatic} 2>> {log}
        """

rule Somaticsniper:
    input:
        disease="results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="results/07_picard/marked_duplicates_{names}_blood.bam"
    output:
        "results/09_variant_calling/somaticsniper/{names}.vcf"
    log:
        "logs/somaticsniper_{names}.log"
    params: 
        ref="-f data/ref_data/hg38.fasta",
        format="-F vcf"
    shell:
        """
        bam-somaticsniper {params.format} {params.ref} {input.disease} {input.blood} {output} 2>> {log}
        """

rule Muse_call:
    input: 
        disease="results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="results/07_picard/marked_duplicates_{names}_blood.bam"
    output: 
        "results/09_variant_calling/MuSE/{names}.MuSE.txt"
    log: 
        "logs/muse_call_{names}.log"
    params:
        ref="-f data/ref_data/hg38.fasta",
        core=" -n 32"
    shell: 
        """
        MuSE call {params.ref} -O results/09_variant_calling/MuSE/{wildcards.names} {params.core} {input.disease} {input.blood} 2>> {log}
        """

rule Muse_sump: 
    input:
        "results/09_variant_calling/MuSE/{names}.MuSE.txt"
    output:
        vcf="results/09_variant_calling/MuSE/{names}.vcf"
    log: 
        "logs/muse_sump_{names}.log"
    params:
        seq="-E",
        core="-n 2"
    shell:
        """
        MuSE sump -I {input} -O {output.vcf} {params.seq} {params.core} 2>> {log}
        """

rule Vardict:
    input: 
        disease="results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="results/07_picard/marked_duplicates_{names}_blood.bam"
    output: 
        "results/09_variant_calling/vardict/{names}.vcf"
    log: 
        "logs/vardict_{names}.log"
    params: 
        ref="-G data/ref_data/hg38.fasta"
    shell: 
        """
        AF_THR="0.01" # minimum allele frequency
        vardict {params.ref} -f $AF_THR -N {wildcards.names} -b "{input.disease}|{input.blood}" | testsomatic.R | var2vcf_paired.pl -N "{wildcards.names}_tumor|{wildcards.names}_blood" -f $AF_THR
        """



