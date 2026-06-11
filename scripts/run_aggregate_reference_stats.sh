#!/bin/bash

# Thin runner for the reference-stats AGGREGATION step (the back half of the
# reference-statistics pipeline). Reads the per-metric reference CSVs computed on
# the MARTS-DB known-TPS set and writes the single committable reference-stats
# JSON. This is the lightweight, CPU-only, no-GPU step — safe to run on a login
# node or laptop; it does NOT run any metric tools (see compute_reference_stats.sh
# for the metric-computation half).

USAGE="--input_dir <dir_of_per_metric_csvs> [--output <json_path>] [--reference_name <name>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --input_dir       Directory holding the per-metric reference CSVs (required)"
    echo "  --output          Output JSON path (optional; default"
    echo "                    src/reference_stats/marts_db_metric_stats.json)"
    echo "  --reference_name  Reference-set label embedded in the JSON (default marts_db)"
    echo "  -h, --help        Show this help message and exit"
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input_dir)
            input_dir="$2"; shift 2 ;;
        --output)
            output="$2"; shift 2 ;;
        --reference_name)
            reference_name="$2"; shift 2 ;;
        -h|--help)
            Help; exit 0 ;;
        *)
            echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$input_dir" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/reference_stats

args=("$input_dir")
[[ -n "$output" ]] && args+=(--output "$output")
[[ -n "$reference_name" ]] && args+=(--reference_name "$reference_name")

python aggregate_reference_stats.py "${args[@]}"
