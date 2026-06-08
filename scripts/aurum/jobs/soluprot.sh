#!/bin/bash
#SBATCH -J SoluProt
#SBATCH --time=0-03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
# Aurum3's custom submit plugin auto-selects the partition; the old `a36_any`
# partition no longer exists. Request the gen-a (a36) hardware via constraint.
#SBATCH --constraint=gen-a

# Usage: sbatch soluprot.sh --fasta_path <fasta_path>

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_soluprot.sh $@
