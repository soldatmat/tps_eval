#!/bin/bash
#SBATCH -J embESM1b
#SBATCH --constraint=gen-b
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch esm_embedding.sh --fasta_path <fasta_path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_esm_embedding.sh $@
