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
    if sample.endswith("_1.fq.gz"):
        samples = sample.split("_1.fq")[0]
        sample_names.append(samples)
        #print(sample_names)

rule all:
    input:

rule fastqc: 
    input:
        "data//{names}_{con}.fq.gz"
    output:
        result = directory("results/00_fastqc/{names}_{con}_fastqc/")
    log:
        "logs/fastqc_{names}_{con}.log"
    conda:
        "envs/fastqc.yaml"
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
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc results/00_fastqc/ -o {output.result} 2>> {log}
        """

rule Fastp:
    input:
        first = "data/{names}_1.fq.gz",
        second = "data/{names}_2.fq.gz"
    output:
        first = "results/04_fastp/{names}_1.fq.gz",
        second = "results/04_fastp/{names}_2.fq.gz",
        html = "results/04_fastp/{names}_fastp.html",
        json = "results/04_fastp/{names}_fastp.json"
    params:
        extra="-w 16"
    log:
        "logs/fastp_{names}.log"
    conda:
        "envs/fastp.yaml"
    shell:
        """
        fastp {params.extra} -i {input.first} -I {input.second} -o {output.first} -O {output.second} -h {output.html} -j {output.json} --detect_adapter_for_pe 2>> {log}
        """