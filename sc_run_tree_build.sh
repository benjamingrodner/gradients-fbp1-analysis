#!/bin/bash

# --- SLURM Settings ---
#SBATCH --job-name=snakemake    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per file
#SBATCH --mem=500G                      # Memory limit
#SBATCH --time=2-00:00:00               # Time limit (hrs:min:sec)
#SBATCH --output=logs/snakemake.log           # Standard output and error log (%j = JobID)

CONDA_PATH=$(conda info --base)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate snakemake_host

snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --jobs 16 \
    --use-conda \
    --printshellcmds \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --resources cores=16 \
    --executor slurm \
    --default-resources \
        slurm_account="bgrodner" \
        slurm_partition="main" \
        mem_mb=$((1 * 1024)) \
        runtime=$((1 * 3600)) \
    --
    
    # -n \