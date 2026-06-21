#!/bin/bash
#PBS -N EnzymeExplorer
#PBS -l walltime=02:30:00
#PBS -l select=1:ncpus=20:ngpus=1:mem=20gb:gpu_mem=20gb:scratch_local=10gb

# Usage: qsub -v args="[--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>" enzyme_explorer.sh

module add mambaforge # run_enzyme_explorer.sh activates the conda environment set in paths.sh


# run_eval_pipeline.py passes the runner argv base64-encoded as $args_b64 (so commas in
# args survive PBS -v); decode to $args. $args (plain) is the manual-submission fallback.
[ -n "$args_b64" ] && args="$(printf %s "$args_b64" | base64 -d)"
# $tps_eval_root is passed by run_eval_pipeline.py via `qsub -v`; fall back to
# $PBS_O_WORKDIR for manual submit_job.sh submission from the repo root.
TPS_EVAL_ROOT="${tps_eval_root:-$PBS_O_WORKDIR}"
. "$TPS_EVAL_ROOT/paths.sh" # load TPS_EVAL_ROOT, ENZYME_EXPLORER_SEQUENCE_ONLY_PATH variables


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export TORCH_HOME="$ENZYME_EXPLORER_SEQUENCE_ONLY_PATH/data/torch_cache"


cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_enzyme_explorer.sh with args: $args"
sh run_enzyme_explorer.sh $args
