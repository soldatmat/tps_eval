#!/bin/bash
#PBS -N EnzymeExplorer
#PBS -l walltime=02:30:00
#PBS -l select=1:ncpus=20:ngpus=1:mem=20gb:gpu_mem=20gb:scratch_local=10gb

# Usage: sbatch enzyme_explorer.sh [--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>

module add mambaforge # run_enzyme_explorer.sh activates the conda environment set in paths.sh

# test if the scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export TORCH_HOME=/storage/brno2/home/soldatmat/documents/terpene_synthases/EnzymeExplorer/data/torch_cache

echo $(pwd)
cd $(dirname "$BASH_SOURCE")/../..
echo $(pwd)

sh run_enzyme_explorer.sh $@
