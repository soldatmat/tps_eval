#!/bin/bash
#SBATCH -J min_embedding_distance
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --partition=a36_96

# Usage: sbatch min_embedding_distance.sh --embeddings_path <embeddings_path> [--train_embeddings_path <train_embeddings_path>]

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_min_embedding_distance.sh $@
