#!/bin/bash
#PBS -N embESM1b
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:ngpus=1:mem=32gb:gpu_mem=16gb:scratch_local=20gb:gpu_cap=compute_70

# Usage: qsub -v "args_b64=<base64 runner args>,tps_eval_root=<repo>" esm_embedding.sh
#   args_b64 is the base64 of the runner argv (run_eval_pipeline.py sets it so
#   commas inside args survive PBS -v parsing). For manual submission you may
#   instead pass plain args via -v "args=..." (no commas) plus tps_eval_root.
# PBS Pro port of scripts/karolina/jobs/esm_embedding.sh (calls run_esm_embedding.sh).

module add mambaforge  # the runner activates the conda env named in paths.sh

[ -n "$args_b64" ] && args="$(printf %s "$args_b64" | base64 -d)"
TPS_EVAL_ROOT="${tps_eval_root:-$PBS_O_WORKDIR}"
. "$TPS_EVAL_ROOT/paths.sh"  # load env names + external-tool/DB paths

test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR
export HF_HOME="$TPS_EVAL_ROOT/.cache/huggingface"
export TORCH_HOME="$TPS_EVAL_ROOT/.cache/torch"

cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_esm_embedding.sh with args: $args"
sh run_esm_embedding.sh $args

clean_scratch
