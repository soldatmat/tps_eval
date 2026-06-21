#!/bin/bash
#PBS -N esmfold
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:ngpus=1:mem=48gb:gpu_mem=24gb:scratch_local=40gb:gpu_cap=compute_70

# Usage: qsub -v "args_b64=<base64 runner args>,tps_eval_root=<repo>" esmfold.sh
#   args_b64 is the base64 of the runner argv (run_eval_pipeline.py sets it so
#   commas inside args survive PBS -v parsing). For manual submission you may
#   instead pass plain args via -v "args=..." (no commas) plus tps_eval_root.
# PBS Pro port of scripts/karolina/jobs/esmfold.sh (calls run_esmfold.sh).

module add mambaforge  # the runner activates the conda env named in paths.sh

[ -n "$args_b64" ] && args="$(printf %s "$args_b64" | base64 -d)"
TPS_EVAL_ROOT="${tps_eval_root:-$PBS_O_WORKDIR}"
. "$TPS_EVAL_ROOT/paths.sh"  # load env names + external-tool/DB paths

test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export HF_HOME="$TPS_EVAL_ROOT/.cache/huggingface"
export TORCH_HOME="$TPS_EVAL_ROOT/.cache/torch"

cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_esmfold.sh with args: $args"
sh run_esmfold.sh $args

clean_scratch
