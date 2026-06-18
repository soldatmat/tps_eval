#!/bin/bash
#SBATCH -J eval_pipeline_continuation
#SBATCH --constraint=gen-a
#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# mg_ee auto-chain continuation: re-invoke the pipeline after EE completes. Light (only
# submits jobs). Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_eval_pipeline_continuation.sh "$@"
