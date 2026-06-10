#!/bin/bash
#SBATCH -J motif_pair_distance
#SBATCH --constraint=gen-a
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G

# Usage: sbatch motif_pair_distance.sh --fasta_path <fasta_path>
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_motif_pair_distance.sh "$@"
