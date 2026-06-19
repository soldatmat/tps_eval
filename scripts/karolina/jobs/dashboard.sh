#!/bin/bash
#SBATCH -J dashboard
#SBATCH --time=0-00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=qcpu

# Builds the interactive natural-bands HTML dashboard (pure Python stdlib; no conda env).
# Usage: sbatch dashboard.sh --designs <glob|dir|csv ...> [--output <out.html>] [--bands <json ...>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_build_dashboard.sh "$@"
