#!/bin/bash
#PBS -N EnzymeExplorer
#PBS -l walltime=02:30:00
#PBS -l select=1:ncpus=20:ngpus=1:mem=20gb:gpu_mem=20gb:scratch_local=10gb

# Usage: qsub -v args="[--sequences_csv_path=<sequences_csv_path> --fasta_path=<fasta_path>] --structs_dir=<structs_dir>" enzyme_explorer.sh

module add mambaforge # run_enzyme_explorer.sh activates the conda environment set in paths.sh


# TODO remove absolute path
. /storage/brno2/home/soldatmat/documents/terpene_synthases/tps_eval/paths.sh # load TPS_EVAL_ROOT, ENZYME_EXPLORER_PATH variables


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export TORCH_HOME="$ENZYME_EXPLORER_PATH/data/torch_cache"


cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_enzyme_explorer.sh with args: $args"
sh run_enzyme_explorer.sh $args
