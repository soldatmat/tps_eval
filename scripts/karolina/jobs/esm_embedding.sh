#!/bin/bash
#SBATCH -J embESM1b
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch esm_embedding.sh --fasta_path <fasta_path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

# Redirect torch's cache off $HOME (tiny ~24 GB quota, usually full) to the shared
# project cache, where the esm1b weights are already cached. Without this, esm
# tries to download 7.8 GB of weights into a full $HOME and fails.
export TORCH_HOME=/mnt/proj2/fta-26-15/.cache/torch

sh run_esm_embedding.sh "$@"
