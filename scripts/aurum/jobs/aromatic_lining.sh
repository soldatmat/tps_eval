#!/bin/bash
#SBATCH -J aromatic_lining
#SBATCH --constraint=gen-a
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Usage: sbatch aromatic_lining.sh --structs_dir <structs_dir> [--save_path <save_path>] [--cutoff <A>] [--cation_pi_min <A>] [--cation_pi_max <A>] [--face_angle_deg <deg>]
# Aurum3's submit plugin auto-selects the partition from --constraint/--time/--mem.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_aromatic_lining.sh "$@"
