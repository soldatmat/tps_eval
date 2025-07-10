#!/bin/bash
#SBATCH -J max_sequence_identity
#SBATCH --time=0-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=a36_any

# Usage: sbatch max_sequence_identity.sh --fasta_path <fasta_path> [--train_path <train_path> --train]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

JULIA_NUM_THREADS=10 sh run_max_sequence_identity.sh "$@"
