#!/bin/bash
#SBATCH -J swissprot_search
#SBATCH --constraint=gen-a
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

# Usage: sbatch swissprot_search.sh --fasta_path <fasta> [--save_path <csv> --top_n <n> --threads <n>]
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.
# DIAMOND blastp parallelizes well; run_swissprot_search.sh defaults --threads to
# $SLURM_CPUS_PER_TASK.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_swissprot_search.sh "$@"
