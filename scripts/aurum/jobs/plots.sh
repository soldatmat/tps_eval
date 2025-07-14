#!/bin/bash
#SBATCH -J plots
#SBATCH --time=0-00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --partition=a36_any

# Usage: sbatch plots.sh --fasta_paths <fasta1.fa> [<fasta2.fa> ...] --data_names <name1> [<name2> ...] --data_colors <color1> [<color2> ...] [--targets <target1> <target2> ... --save_dir <save_dir>]"


SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')

cd $(dirname "$SCRIPT_PATH")/../..

sh run_plots.sh $@
