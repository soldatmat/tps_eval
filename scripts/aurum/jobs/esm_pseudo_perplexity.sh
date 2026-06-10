#!/bin/bash
#SBATCH -J esmppl
#SBATCH --constraint=gen-b
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch esm_pseudo_perplexity.sh --fasta_path <fasta_path> [--method swoop|masked]
# Light naturalness metric (ESM-1b masked-LM). Use gen-b RTX 3090, matching the other
# ESM/MPNN GPU jobs: gen-a + gpu:1 routes to the single-node a36_96_gpu partition
# (node a233), which is frequently down and leaves the job PENDING (and blocks plots).
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem/--gres.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_esm_pseudo_perplexity.sh "$@"
