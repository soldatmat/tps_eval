#!/bin/bash
#SBATCH -J aggregation
#SBATCH --time=0-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=qcpu

# Usage: sbatch aggregation.sh --structs_dir <structs_dir> [--save_path <save_path>] [--save_residue_scores] [--residue_scores_dir <dir>]
# Aggrescan3D static mode is CPU-only (one structure at a time, ~seconds each).

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_aggregation.sh "$@"
