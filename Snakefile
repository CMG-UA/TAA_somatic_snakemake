import os
import csv

CONDITIONS = ["1", "2"]
TYPE = ["tumor", "blood"]

# Define directories
#REFDIR = os.getcwd()
#print(REFDIR)
#sample_dir = REFDIR+"/data/"
sample_dir = "/home/mhannaert/TAA_somatic_snakemake/data/"

sample_names = []
sample_list = os.listdir(sample_dir)
for i in range(len(sample_list)):
    sample = sample_list[i]
    if sample.endswith("_tumor_1.fastq.gz"):
        samples = sample.split("_tumor_1.fastq")[0]
        sample_names.append(samples)
        #print(sample_names)

known_genes_list = []
with open("/home/mhannaert/TAA_somatic_snakemake/data/known_genes.txt", "r") as known_genes:
    reader = csv.reader(known_genes)
    next(reader)  # sla de header over: "symbol,chromosoom,start,stop"
    for row in reader:
        symbol, chrom, start, stop = row
        known_genes_list.append({
            "symbol": symbol,
            "chromosome": chrom,
            "start": int(start),
            "stop": int(stop)
        })


hematopoëse_genes_list = []
with open("/home/mhannaert/TAA_somatic_snakemake/data/known_genes.txt", "r") as hema_genes:
    reader = csv.reader(hema_genes)
    next(reader)  # sla header over
    for row in reader:
        symbol, chrom, start, stop = row
        hematopoëse_genes_list.append({
            "symbol": symbol,
            "chromosome": chrom,
            "start": int(start),
            "stop": int(stop)
        })

rule all:
    input:"/home/mhannaert/TAA_somatic_snakemake/results/01_multiqc/multiqc_report.html",
        "/home/mhannaert/TAA_somatic_snakemake/results/04_ATmultiqc/multiqc_report.html",
       #expand("results/09_variant_calling/mutect/{names}_somatic.vcf.gz", names=sample_names),
        #expand("results/09_variant_calling/varscan/{names}.vcf.gz", names=sample_names),
        #expand("results/09_variant_calling/somaticsniper/{names}.vcf.gz", names=sample_names),
        expand("/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/germline/{names}_selected_germline.vcf", names=sample_names),
        expand("/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/{names}_consensus.vcf", names=sample_names)
        #expand("results/09_variant_calling/MuSE/{names}.vcf", names=sample_names)
        

rule fastqc: 
    input:
        "/home/mhannaert/TAA_somatic_snakemake/data/{names}_{type}_{con}.fastq.gz"
    output:
        result = directory("/home/mhannaert/TAA_somatic_snakemake/results/00_fastqc/{names}_{type}_{con}_fastqc/")
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/fastqc_{names}_{type}_{con}.log"
    params:
        extra="-t 64"
    shell:
        """
        fastqc {params.extra} {input} --extract -o /home/mhannaert/TAA_somatic_snakemake/results/00_fastqc/ 2>> {log}
        """

rule multiqc:
    input:
        expand("/home/mhannaert/TAA_somatic_snakemake/results/00_fastqc/{names}_{type}_{con}_fastqc/", names=sample_names, type=TYPE, con = CONDITIONS)
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/01_multiqc/multiqc_report.html",
        result = directory("/home/mhannaert/TAA_somatic_snakemake/results/01_multiqc/")  
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/multiqc.log"
    shell:
        """
        multiqc /home/mhannaert/TAA_somatic_snakemake/results/00_fastqc/ -o {output.result} 2>> {log}
        """

rule Fastp:
    input:
        first = "/home/mhannaert/TAA_somatic_snakemake/data/{names}_{type}_1.fastq.gz",
        second = "/home/mhannaert/TAA_somatic_snakemake/data/{names}_{type}_2.fastq.gz"
    output:
        first = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_1.fastq.gz",
        second = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_2.fastq.gz",
        html = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_fastp.html",
        json = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_fastp.json"
    params:
        extra="-w 16"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/fastp_{names}_{type}.log"
    shell:
        """
        fastp {params.extra} -i {input.first} -I {input.second} -o {output.first} -O {output.second} -h {output.html} -j {output.json} --detect_adapter_for_pe 2>> {log}
        """
rule aftertrimming_fastqc: 
    input:
        "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_{con}.fastq.gz"
    output:
        result = directory("/home/mhannaert/TAA_somatic_snakemake/results/03_ATfastqc/{names}_{type}_{con}_fastqc/")
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/ATfastqc_{names}_{type}_{con}.log"
    params:
        extra="-t 64"
    shell:
        """
        fastqc {params.extra} {input} --extract -o /home/mhannaert/TAA_somatic_snakemake/results/03_ATfastqc/ 2>> {log}
        """

rule aftertrimming_multiqc:
    input:
        expand("/home/mhannaert/TAA_somatic_snakemake/results/03_ATfastqc/{names}_{type}_{con}_fastqc/", names=sample_names, type=TYPE, con = CONDITIONS)
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/04_ATmultiqc/multiqc_report.html",
        result = directory("/home/mhannaert/TAA_somatic_snakemake/results/04_ATmultiqc/")  
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/ATmultiqc.log"
    shell:
        """
        multiqc /home/mhannaert/TAA_somatic_snakemake/results/03_ATfastqc/ -o {output.result} 2>> {log}
        """

rule BWA_indexing: 
    input: 
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta"
    output: 
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.amb",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.ann",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.bwt",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.pac",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.sa"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/bwa_index.log"
    shell: 
        """
        bwa index {input} 2>> {log}
        """

rule BWA_alignment:
    input: 
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.amb",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.ann",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.bwt",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.pac",
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta.sa",
        first = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_1.fastq.gz",
        second = "/home/mhannaert/TAA_somatic_snakemake/results/02_fastp/{names}_{type}_2.fastq.gz"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/05_bwa/{names}_{type}.sam"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/bwa_{names}_{type}.log"
    params:
        extra="-t 64"
    shell:
     """
     bwa mem {params.extra} /home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta {input.first} {input.second} > {output} 2>> {log}
     """

rule samtools: 
    input:
        "/home/mhannaert/TAA_somatic_snakemake/results/05_bwa/{names}_{type}.sam"
    output: 
        "/home/mhannaert/TAA_somatic_snakemake/results/06_samtools/{names}_{type}.bam"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/samtools_{names}_{type}.log"
    params:
        threads="-@ 64"
    shell:
        """
        samtools sort {params.threads} {input} -o {output} 2>> {log}
        samtools index {params.threads} {output} 2>> {log}
        """ 

rule AddReadGroups:
    input:
        "/home/mhannaert/TAA_somatic_snakemake/results/06_samtools/{names}_{type}.bam"
    output:
        first="/home/mhannaert/TAA_somatic_snakemake/results/06_samtools/{names}_{type}_rg.bam",
        second="/home/mhannaert/TAA_somatic_snakemake/results/06_samtools/{names}_{type}_rg_sorted.bam"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/add_read_groups_{names}_{type}.log"
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
        "/home/mhannaert/TAA_somatic_snakemake/results/06_samtools/{names}_{type}_rg_sorted.bam"
    output:
        bam="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_{type}.bam",
        txt="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_dup_metrics_{names}_{type}.txt"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/picard_{names}_{type}.log"
    params:
        regex="--READ_NAME_REGEX null",
        resources="--MAX_RECORDS_IN_RAM 2000000 --COMPRESSION_LEVEL 5"
    shell: 
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.txt} {params.regex} {params.resources} 2>> {log}
        """
rule mpileup:
    input: 
        "/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_{type}.bam"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/08_mpileup/marked_duplicates_{names}_{type}.pileup"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/mpileup_{names}_{type}.log"
    params: 
        ref="-f /home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta"
    shell: 
        """
        samtools mpileup {params.ref} {input} > {output} 2>> {log}
        """

rule Varscan: 
    input: 
        disease="/home/mhannaert/TAA_somatic_snakemake/results/08_mpileup/marked_duplicates_{names}_tumor.pileup",
        blood="/home/mhannaert/TAA_somatic_snakemake/results/08_mpileup/marked_duplicates_{names}_blood.pileup"
    output: 
        snp="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.snp",
        indel="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.indel"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/varscan_{names}.log"
    #params: 
    shell: 
        """
        varscan somatic {input.blood} {input.disease} --output-snp {output.snp} --output-indel {output.indel} 2>> {log}
        """

rule snp_to_vcf:
    input:
        snp="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.snp"
    output:
        first="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf",
        second="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/convert_snp_vcf_{names}.log"
    shell:
        """
        python3 /home/mhannaert/TAA_somatic_snakemake/scripts/vs_format_converter.py {input.snp} >> {output.first} 2>> {log}
        bgzip -k {output.first}
        """
rule Getting_germline_variants: 
    input:
        vcf="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/germline/{names}_germline.vcf.gz"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/getting_germline_variants_{names}.log"
    shell:
        """
        bcftools view -i 'SS="1"' {input} -Oz -o {output} 2>> {log}
        """
rule selecting_germline_variants:
    input:
        germ_file="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/germline/{names}_germline.vcf.gz",
        bed="/home/mhannaert/TAA_somatic_snakemake/data/known_genes.bed"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/germline/{names}_selected_germline.vcf"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/selecting_germline_variants_{names}.log"
    shell:
        """
        zcat {input.germ_file} | bedtools intersect -a - -b {input.bed} -header > {output} 2>> {log}
        """
rule fasta_dict: 
    input:
        ref="/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.dict"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/fasta_dict.log"
    shell:
        """
        gatk CreateSequenceDictionary -R {input.ref} 2>> {log}
        """
        
rule index:
    input:
        "/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_{type}.bam"
    output: 
        "/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_{type}.bam.bai"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/index_dup_{names}_{type}.log"
    params: 
        "-@ 32"
    shell: 
        "samtools index {params} {input} 2>> {log}"

rule Mutect: 
    input: 
        disease="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_blood.bam",
        index_d="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_tumor.bam.bai",
        index_b="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_blood.bam.bai",
        ref_dict="/home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.dict"
    output: 
       "/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/mutect/{names}_somatic.vcf.gz"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/mutect_{names}.log"
    params: 
        ref="-R /home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta",
        threads="--native-pair-hmm-threads 16"
    shell: 
        """
        gatk Mutect2 \
        {params.threads} \
        {params.ref} \
        -I {input.disease}\
        -I {input.blood} \
        -normal {wildcards.names}_blood \
        -O {output} 2>> {log}
        """

rule Somaticsniper:
    input:
        disease="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_tumor.bam",
        blood="/home/mhannaert/TAA_somatic_snakemake/results/07_picard/marked_duplicates_{names}_blood.bam"
    output:
        first="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf",
        second="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf.gz"
    log:
        "/home/mhannaert/TAA_somatic_snakemake/logs/somaticsniper_{names}.log"
    params: 
        ref="-f /home/mhannaert/TAA_somatic_snakemake/data/ref_data/hg38.fasta",
        format="-F vcf",
        quality="-Q 40 -G -L" #recommended quality parameters for mapping with bwa
    shell:
        """
        bam-somaticsniper {params.quality} {params.format} {params.ref} {input.disease} {input.blood} {output.first} 2>> {log}
        bgzip -k {output.first} 2>> {log}
        """

rule index_tbi:
    input:
        mutect="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/mutect/{names}_somatic.vcf.gz",
        somaticsniper="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf.gz",
        varscan="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz"
    output:
        mutect="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/mutect/{names}_somatic.vcf.gz.tbi",
        somaticsniper="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf.gz.tbi",
        varscan="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz.tbi"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/index_tbi_{names}.log"
    params:
        threads="-@ 32",
        filetype="-p vcf"
    shell:
        """
        tabix {params.threads} {params.filetype} {input.mutect} &
        tabix {params.threads} {params.filetype} {input.somaticsniper} &
        tabix {params.threads} {params.filetype} {input.varscan} 2>> {log}
        """

rule consensus_bcftools:
    input:
        mutect="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/mutect/{names}_somatic.vcf.gz",
        somaticsniper="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf.gz",
        varscan="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz",
        mutect_tbi="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/mutect/{names}_somatic.vcf.gz.tbi",
        somaticsniper_tbi="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/somaticsniper/{names}.vcf.gz.tbi",
        varscan_tbi="/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/varscan/{names}.vcf.gz.tbi"
    output:
        "/home/mhannaert/TAA_somatic_snakemake/results/09_variant_calling/{names}_consensus.vcf"
    log: 
        "/home/mhannaert/TAA_somatic_snakemake/logs/consensus_bcf_{names}.log"
    params: 
        appear_number= "-n+2", #the -n+2 is for variants in at least 2 out of 3
        metadata= "-w1", #without metadata, only variants records
        threads="--threads 32"
    shell: 
        """
        bcftools isec {params.threads} {params.appear_number} {params.metadata} -o {output} {input.mutect} {input.somaticsniper} {input.varscan} 2>> {log}
        """
