#!/bin/bash
#SBATCH -J TmProt
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch tmprot.sh --fasta_path <fasta_path> [--device cpu|cuda]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

# Redirect the HF/torch caches off $HOME (tiny quota, usually full) to the shared
# project cache, where the esm2_t33_650M_UR50D weights are already cached. Without
# this, transformers tries to download ~2.5 GB of weights into a full $HOME and fails.
export HF_HOME=/mnt/proj2/fta-26-15/.cache/huggingface
export TORCH_HOME=/mnt/proj2/fta-26-15/.cache/torch

sh run_tmprot.sh "$@"
