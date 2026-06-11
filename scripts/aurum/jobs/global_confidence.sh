#!/bin/bash
#SBATCH -J global_confidence
#SBATCH --constraint=gen-a
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Usage: sbatch global_confidence.sh --pae_dir <dir> [--structs_dir <dir>] [--save_path <csv>]
# Global fold confidence (pTM/ipTM) read from the fold-time <ID>_pae.npz files
# (numpy only; CPU). Aurum3's submit plugin auto-selects the partition from
# --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_global_confidence.sh "$@"
