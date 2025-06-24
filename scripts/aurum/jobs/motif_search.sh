#!/bin/bash
#SBATCH -J motif_search
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --partition=a36_96

# Usage: sbatch motif_search.sh --fasta_path <fasta_path> [<motif1> <motif2> ...]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_motif_search.sh $@
