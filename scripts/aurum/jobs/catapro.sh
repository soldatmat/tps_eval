#!/bin/bash
#SBATCH -J CataPro
#SBATCH --constraint=gen-b
#SBATCH --time=0-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:geforce_rtx_3090:1

# Usage: sbatch catapro.sh --fasta_path <fasta_path> [--target_substrate <CODE>] [--out_suffix <s>]
# CataPro kinetics (kcat/Km/kcat-over-Km) via ProtT5+MolT5 -> needs a GPU. Use gen-b
# RTX 3090, matching the other ESM/MPNN GPU jobs: gen-a + gpu:1 routes to the
# single-node a36_96_gpu partition (node a233), which is frequently down and leaves
# the job PENDING. Aurum3's submit plugin auto-selects the partition from
# --constraint/--time/--mem/--gres (do NOT pass -p).

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")/../..

sh run_catapro.sh "$@"
