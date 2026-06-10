#!/bin/bash
#SBATCH -J esmfold
#SBATCH --time=0-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch esmfold.sh --fasta_path <fasta_path> --save_dir <save_dir> [--chunk_size <n>]
# ESMFold is Karolina's folding option (AlphaFold is Aurum-only).

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd "$(dirname "$SCRIPT_PATH")/../.."

# Redirect HuggingFace's cache off $HOME (tiny ~24 GB quota, usually full) to the
# shared project cache. The ESMFold weights (~3 GB) download here once; without
# this, transformers tries to write into a full $HOME and fails.
export HF_HOME=/mnt/proj2/fta-26-15/.cache/huggingface

sh run_esmfold.sh "$@"
