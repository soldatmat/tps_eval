#!/bin/bash
#SBATCH -J knn_label_transfer
#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --constraint=gen-a

# Usage: sbatch knn_label_transfer.sh <calibrate|predict> [--sequence_topk ...] [--embedding_topk ...] [--structural_topk ...] --label_file <csv> --out <path> [--calibration <json>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_knn_label_transfer.sh "$@"
