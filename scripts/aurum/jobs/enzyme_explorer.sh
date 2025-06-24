#!/bin/bash
#SBATCH -J EnzymeExplorer
#SBATCH --partition=b32_128_gpu
#SBATCH --time=03:59:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:2

# Usage: sbatch enzyme_explorer.sh [--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>

# original partition b32_128_gpu
# time = 1:59:00 for (?)1600 sequences

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_enzyme_explorer.sh $@
