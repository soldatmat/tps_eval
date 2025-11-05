#!/bin/bash
#PBS -N EnzymeExplorer_sequence_only
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:ngpus=1:mem=20gb:gpu_mem=20gb:scratch_local=10gb

# Usage: qsub -v args="--fasta_path <fasta_path>" enzyme_explorer_sequence_only.sh

module add mambaforge # run_enzyme_explorer_sequence_only.sh activates the conda environment set in paths.sh


# TODO remove absolute path
. /storage/brno2/home/soldatmat/documents/terpene_synthases/tps_eval/paths.sh # load TPS_EVAL_ROOT, ENZYME_EXPLORER_PATH variables


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export TORCH_HOME="$ENZYME_EXPLORER_PATH/data/torch_cache"


cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_enzyme_explorer_sequence_only.sh with args: $args"
sh run_enzyme_explorer_sequence_only.sh $args
