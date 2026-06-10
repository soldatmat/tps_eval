#!/bin/bash
#SBATCH -J esmppl
#SBATCH --constraint=gen-a
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --gres=gpu:1

# Usage: sbatch esm_pseudo_perplexity.sh --fasta_path <fasta_path> [--method swoop|masked]
# Light naturalness metric (ESM-1b masked-LM): gen-a is fine. Aurum3's submit
# plugin auto-selects the partition from --constraint/--time/--mem/--gres.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_esm_pseudo_perplexity.sh "$@"
