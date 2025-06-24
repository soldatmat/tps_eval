#!/bin/bash
#SBATCH -J EnzymeExplorer_sequence_only
#SBATCH --time=0-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --gres=gpu:1
#SBATCH --partition=b32_128_gpu

# Usage: sbatch enzyme_explorer_sequence_only.sh --fasta_path <fasta_path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_enzyme_explorer_sequence_only.sh $@
