#!/bin/bash
#SBATCH -J domain_structural_identity
#SBATCH --time=0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=qcpu

# Usage: sbatch domain_structural_identity.sh --structs_dir <dir> [--known_domain_structures_root <dir>] [--save_path <csv>] [--n_jobs <n>] [--n_iters <n>]
# EE TPS-domain detection (PyMOL + foldseek) + foldseek domain alignment — CPU-only; no GPU.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_domain_structural_identity.sh "$@"
