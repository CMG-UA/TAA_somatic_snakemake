#!/bin/bash
#SBATCH --job-name=Snakemake_TAA_somatic
#SBATCH --mem=20GB
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=Snakemake_TAA_somatic_%j.out
#SBATCH --error=Snakemake_TAA_somatic_%j.err

#conda activate TAA_somatic_pipeline
source $HOME/.bashrc

/opt/software/miniconda3/bin/snakemake --snakefile /home/mhannaert/TAA_somatic_snakemake/Snakefile --use-conda --cluster "sbatch --mem=40G -t 800 --cpus-per-task=4 --output=log_sbatch/$SLURM_JOB_ID.%j.out --error=log_sbatch/$SLURM_JOB_ID.%j.err" -j16 --verbose --rerun-incomplete --latency-wait 36000

# You can also use env var set by SLURM, such as $SLURM_JOB_ID, $SLURM_NODELIST, etc.