#!/bin/bash
#SBATCH -J EnzymeExplorer
#SBATCH --time=0-02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch enzyme_explorer.sh [--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_enzyme_explorer.sh "$@"
