#!/bin/bash
#SBATCH -J SoluProt
#SBATCH --time=0-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --partition=a36_any

# Usage: sbatch soluprot.sh --fasta_path <fasta_path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_soluprot.sh $@
