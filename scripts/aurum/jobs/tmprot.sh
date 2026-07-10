#!/bin/bash
#SBATCH -J TmProt
#SBATCH --constraint=gen-b
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch tmprot.sh --fasta_path <fasta_path> [--device cpu|cuda]
# TmProt (ESM2-650M + LoRA) melting-temperature predictor. Use gen-b RTX 3090,
# matching the other ESM GPU jobs: gen-a + gpu:1 routes to the single-node
# a36_96_gpu partition (node a233), which is frequently down and leaves the job
# PENDING. Aurum3's submit plugin auto-selects the partition from
# --constraint/--time/--mem/--gres (do NOT pass -p).

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh tool_wrappers/run_tmprot.sh "$@"
