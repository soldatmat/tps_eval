#!/bin/bash
#SBATCH -J mpnnscore
#SBATCH --constraint=gen-b
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch proteinmpnn_score.sh --structs_dir <structs_dir> [--save_path <save_path>]
# ProteinMPNN NLL of each design's own sequence given its backbone. Light model;
# gen-b RTX 3090 schedules fast. Aurum3's submit plugin picks the partition from
# --constraint/--time/--mem/--gres.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_proteinmpnn_score.sh "$@"
