#!/bin/bash
#SBATCH -J pocket_descriptors
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=qcpu

# Usage: sbatch pocket_descriptors.sh --structs_dir <structs_dir> [--save_path <save_path>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_pocket_descriptors.sh "$@"
