#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>] [--save_residue_scores] [--residue_scores_dir <dir>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir           Directory of structures (.pdb/.cif); file stem = ID (required)"
    echo "  --save_path             Output CSV path (optional; default <structs_dir>_aggregation.csv)"
    echo "  --save_residue_scores   Also dump per-residue A3D scores to a side dir (optional; off by default)"
    echo "  --residue_scores_dir    Directory for per-residue scores (optional; default <structs_dir>_aggregation_residue_scores)"
    echo "  -h, --help              Show this help message and exit"
    echo
}

# Parse long options manually
save_residue_scores=0
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir)
            structs_dir="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
            shift 2
            ;;
        --save_residue_scores)
            save_residue_scores=1
            shift 1
            ;;
        --residue_scores_dir)
            residue_scores_dir="$2"
            shift 2
            ;;
        -h|--help)
            Help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            Help
            exit 1
            ;;
    esac
done

if [ -z "$structs_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert structs_dir to absolute path if relative
if [[ "$structs_dir" != /* ]]; then
    structs_dir="$(cd "$structs_dir" && pwd)"
fi
# Convert save_path to absolute path if relative
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load AGGRESCAN3D_ENV

eval "$(conda shell.bash hook)"
conda activate "$AGGRESCAN3D_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/structure_metrics

args=("$structs_dir")
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi
if [[ "$save_residue_scores" == "1" ]]; then
    args+=(--save_residue_scores)
fi
if [[ -n "$residue_scores_dir" ]]; then
    args+=(--residue_scores_dir "$residue_scores_dir")
fi

python run_aggregation.py "${args[@]}"
