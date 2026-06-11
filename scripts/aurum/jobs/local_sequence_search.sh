#!/bin/bash
#SBATCH -J local_sequence_search
#SBATCH --constraint=gen-a
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

# Usage: sbatch local_sequence_search.sh --fasta_path <fasta> [--train_path <ref> --backend <mmseqs2|diamond> --top_k <N>]
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem (no -p).
# Both MMseqs2 and DIAMOND parallelize well; run_local_sequence_search.sh defaults
# --threads to $SLURM_CPUS_PER_TASK.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_local_sequence_search.sh "$@"
