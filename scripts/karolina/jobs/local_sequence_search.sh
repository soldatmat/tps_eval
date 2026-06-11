#!/bin/bash
#SBATCH -J local_sequence_search
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=qcpu

# Usage: sbatch local_sequence_search.sh --fasta_path <fasta> [--train_path <ref> --backend <mmseqs2|diamond> --top_k <N>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_local_sequence_search.sh "$@"
