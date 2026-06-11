#!/bin/bash
#SBATCH -J pocket_descriptors
#SBATCH --constraint=gen-a
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Usage: sbatch pocket_descriptors.sh --structs_dir <structs_dir> [--save_path <save_path>]
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.
# CPU tool (fpocket + P2Rank); --constraint=gen-a (no -p).

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_pocket_descriptors.sh "$@"
