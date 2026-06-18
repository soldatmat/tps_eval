#!/bin/bash
#SBATCH -J eval_pipeline_continuation
#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=qcpu

# mg_ee auto-chain continuation: re-invoke the pipeline after EE completes. Light (only
# submits jobs). NOTE: AlphaFold3 is Aurum-only, so the mg_ee auto-chain (which requires AF3
# folding) is cap-gated off on Karolina; this script exists for symmetry / future ports.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_eval_pipeline_continuation.sh "$@"
