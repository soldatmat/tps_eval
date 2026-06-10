#!/bin/bash
#SBATCH -J domain_composition
#SBATCH --constraint=gen-a
#SBATCH --time=0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G

# Usage: sbatch domain_composition.sh --structs_dir <dir> [--save_path <csv>] [--detections_json <json>] [--n_jobs <n>] [--n_iters <n>]
# EnzymeExplorer TPS-domain detection is PyMOL + foldseek (CPU-only; no GPU).
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_domain_composition.sh "$@"
