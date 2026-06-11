#!/bin/bash
#SBATCH -J extract_pae
#SBATCH --constraint=gen-a
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Usage: sbatch extract_pae.sh --structs_dir <dir> --pae_dir <dir>
# CPU-only (numpy/json). Aurum3's submit plugin auto-selects the partition from
# --constraint/--time/--mem. Post-fold step of the orchestrator's --fold alphafold3 path.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_extract_pae.sh "$@"
