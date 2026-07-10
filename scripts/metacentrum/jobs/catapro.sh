#!/bin/bash
#PBS -N CataPro
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:ngpus=1:mem=32gb:gpu_mem=16gb:scratch_local=20gb:gpu_cap=compute_70

# Usage: qsub -v "args_b64=<base64 runner args>,tps_eval_root=<repo>" catapro.sh
#   args_b64 is the base64 of the runner argv (run_eval_pipeline.py sets it so
#   commas inside args survive PBS -v parsing). For manual submission you may
#   instead pass plain args via -v "args=..." (no commas) plus tps_eval_root.
# PBS Pro port of scripts/karolina/jobs/catapro.sh (calls run_catapro.sh). CataPro
# kinetics via ProtT5+MolT5 -> needs a GPU. The backbones are local
# (vendor/CataPro/models), so no HF/torch weight download is needed.

module add mambaforge  # the runner activates the conda env named in paths.sh

[ -n "$args_b64" ] && args="$(printf %s "$args_b64" | base64 -d)"
TPS_EVAL_ROOT="${tps_eval_root:-$PBS_O_WORKDIR}"
. "$TPS_EVAL_ROOT/paths.sh"  # load env names + external-tool/DB paths

test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
export TMPDIR=$SCRATCHDIR

cd "$TPS_EVAL_ROOT/scripts"
echo "Calling run_catapro.sh with args: $args"
sh tool_wrappers/run_catapro.sh $args

clean_scratch
