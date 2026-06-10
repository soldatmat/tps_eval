#!/bin/bash
#SBATCH -J esmfold
#SBATCH --constraint=gen-b
#SBATCH --time=0-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch esmfold.sh --fasta_path <fasta_path> --save_dir <save_dir> [--chunk_size <n>]
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem/--gres.
# gen-b RTX 3090 (24 GB) preferred over gen-a's single contended Quadro RTX 5000 (16 GB):
# schedules faster and has enough VRAM for long (>500 aa) TPS sequences.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

# Cache the ESMFold weights (~3 GB ESM-2 backbone + folding head) under a project
# dir, NOT $HOME. Edit to a NFS-IB project path; falls back to the default cache.
export HF_HOME="${HF_HOME:-/home/soldat/documents/.cache/huggingface}"

sh run_esmfold.sh "$@"
