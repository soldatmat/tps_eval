#!/bin/bash
#PBS -N knn_label_transfer
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=4gb

# Usage: qsub -v "args_b64=<base64 runner args>,tps_eval_root=<repo>" knn_label_transfer.sh
#   args_b64 is the base64 of the runner argv (run_eval_pipeline.py sets it so
#   commas inside args survive PBS -v parsing). For manual submission you may
#   instead pass plain args via -v "args=..." (no commas) plus tps_eval_root.
# PBS Pro port of scripts/karolina/jobs/knn_label_transfer.sh (calls run_knn_label_transfer.sh).

module add mambaforge  # the runner activates the conda env named in paths.sh

[ -n "$args_b64" ] && args="$(printf %s "$args_b64" | base64 -d)"
TPS_EVAL_ROOT="${tps_eval_root:-$PBS_O_WORKDIR}"
. "$TPS_EVAL_ROOT/paths.sh"  # load env names + external-tool/DB paths

test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR

cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_knn_label_transfer.sh with args: $args"
sh run_knn_label_transfer.sh $args

clean_scratch
