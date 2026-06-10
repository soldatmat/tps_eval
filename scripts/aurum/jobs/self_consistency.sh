#!/bin/bash
#SBATCH -J scrmsd
#SBATCH --constraint=gen-b
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch self_consistency.sh --structs_dir <structs_dir> [--num_seqs <n>] [--limit <n>] [--ids <id...>]
# Self-consistency scRMSD: ProteinMPNN sampling -> ESMFold refold -> Ca-RMSD.
# Heavy (N folds per structure); gen-b RTX 3090 (24 GB) holds ESMFold for TPS-length
# seqs. Aurum3's submit plugin picks the partition from --constraint/--time/--mem/--gres.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

# ESMFold weights cache off $HOME (set HF_HOME, same as esmfold.sh).
export HF_HOME="${HF_HOME:-/home/soldat/documents/.cache/huggingface}"

sh run_self_consistency.sh "$@"
