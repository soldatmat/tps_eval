#!/bin/bash
#SBATCH -J aggregate_reference_stats
#SBATCH --constraint=gen-a
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G

# Usage: sbatch aggregate_reference_stats.sh --input_dir <dir> [--output <json>] [--reference_name <name>]
# Back half of the reference-statistics pipeline: read the per-metric reference
# CSVs (computed on the MARTS-DB known-TPS set) and write the single committable
# reference-stats JSON. CPU-only, light. Aurum3's submit plugin auto-selects the
# partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_aggregate_reference_stats.sh "$@"
