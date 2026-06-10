#!/bin/bash
#SBATCH -J structural_identity
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=qcpu

# Usage: sbatch structural_identity.sh --structs_dir <dir> --known_structs_dir <dir> [--save_path <csv>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_structural_identity.sh "$@"
