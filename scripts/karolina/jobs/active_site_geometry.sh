#!/bin/bash
#SBATCH -J active_site_geometry
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=qcpu

# Usage: sbatch active_site_geometry.sh --structs_dir <structs_dir> [--save_path <save_path>] [--templates <ID,ID,...>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_active_site_geometry.sh "$@"
