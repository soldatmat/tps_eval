#!/bin/bash
#SBATCH -J scrmsd
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch self_consistency.sh --structs_dir <structs_dir> [--num_seqs <n>] [--limit <n>] [--ids <id...>]
# Self-consistency scRMSD: ProteinMPNN sampling -> ESMFold refold -> Ca-RMSD.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

# ESMFold weights cache off $HOME (set HF_HOME, same as Karolina esmfold.sh).
export HF_HOME=/mnt/proj2/fta-26-15/.cache/huggingface

sh run_self_consistency.sh "$@"
