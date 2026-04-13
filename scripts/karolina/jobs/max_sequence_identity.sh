#!/bin/bash
#SBATCH -J max_sequence_identity
#SBATCH --time=0-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=qcpu

# Usage: sbatch max_sequence_identity.sh --fasta_path <fasta_path> [--train_path <train_path> --train]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
sh run_max_sequence_identity.sh "$@"
