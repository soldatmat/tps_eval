#!/bin/bash
#SBATCH -J sdr_divergence
#SBATCH --constraint=gen-a
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Usage: sbatch sdr_divergence.sh --structs_dir <s> --known_structs_dir <k> [--structural_topk <csv>] [--sequence_topk <csv>] [--sdr_panel <csv>] ...
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_sdr_divergence.sh "$@"
