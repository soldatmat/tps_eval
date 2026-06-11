#!/bin/bash
#SBATCH -J substrate_class
#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --constraint=gen-a

# Usage: sbatch substrate_class.sh [--sequence_topk ...] [--embedding_topk ...] [--structural_topk ...] --label_file <csv> --calibration <json> [--pocket_csv <csv>] [--ee_csv <csv>] --out <path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh run_substrate_class.sh "$@"
