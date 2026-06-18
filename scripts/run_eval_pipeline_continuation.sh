#!/bin/bash

# Continuation runner for the `--af3_cofold mg_ee` AUTO-CHAIN. Queued as a job that runs
# AFTER the in-pipeline enzyme_explorer_sequence_only job has produced the EE CSV, then
# re-invokes run_eval_pipeline.py with --enzymeexplorer_csv now present -- so the AF3 fan-out
# co-folds each design with its EE-predicted substrate and submits the structure branch. It
# only SUBMITS jobs (light, login-node-class); no heavy compute.
#
# Args: <orig_cwd> <original run_eval_pipeline.py argv...> --enzymeexplorer_csv <ee_csv>
# We cd to <orig_cwd> first so relative paths in the forwarded argv resolve as they did at
# launch (sbatch must be callable from the compute node — true on Aurum).

if [[ $# -lt 1 ]]; then echo "Usage: $0 <orig_cwd> <run_eval_pipeline argv...>"; exit 1; fi
orig_cwd="$1"; shift

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
REPO=$(pwd)
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"

cd "$orig_cwd"
echo "[continuation] EE complete; re-invoking the pipeline from $orig_cwd (mg_ee inline)"
python "$REPO/scripts/run_eval_pipeline.py" "$@"
