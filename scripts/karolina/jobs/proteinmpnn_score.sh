#!/bin/bash
#SBATCH -J mpnnscore
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch proteinmpnn_score.sh --structs_dir <structs_dir> [--save_path <save_path>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_proteinmpnn_score.sh "$@"
