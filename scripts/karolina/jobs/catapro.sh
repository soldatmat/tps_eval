#!/bin/bash
#SBATCH -J CataPro
#SBATCH --time=0-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=qgpu

# Usage: sbatch catapro.sh --fasta_path <fasta_path> [--target_substrate <CODE>] [--out_suffix <s>]
# CataPro kinetics via ProtT5+MolT5 -> needs a GPU (qgpu). The ProtT5/MolT5 backbones
# are local (vendor/CataPro/models), so no weight download / cache redirect is needed.

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd "$(dirname "$SCRIPT_PATH")/../.."

sh tool_wrappers/run_catapro.sh "$@"
