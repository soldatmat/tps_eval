#!/bin/bash
#SBATCH -J interdomain_pae
#SBATCH --constraint=gen-a
#SBATCH --time=0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G

# Usage: sbatch interdomain_pae.sh --structs_dir <dir> --pae_dir <dir> [--save_path <csv>] [--detections_json <json>] [--per_pair] [--n_jobs <n>] [--n_iters <n>]
# Inter-domain PAE: EE TPS-domain detection (PyMOL + foldseek, CPU-only; no GPU) +
# reduction over the fold-time PAE matrices. Aurum3's submit plugin auto-selects the
# partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_interdomain_pae.sh "$@"
