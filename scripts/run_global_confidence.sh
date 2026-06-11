#!/bin/bash

USAGE="--pae_dir <pae_dir> [--structs_dir <structs_dir>] [--save_path <save_path>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --pae_dir       Directory of <ID>_pae.npz matrices saved at fold time; ID = stem (required)"
    echo "  --structs_dir   Structs dir used only to NAME the output CSV (<structs_dir>_global_confidence.csv) (optional)"
    echo "  --save_path     Output CSV path (optional; default <structs_dir|pae_dir>_global_confidence.csv)"
    echo "  -h, --help      Show this help message and exit"
    echo
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --pae_dir)
            pae_dir="$2"
            shift 2
            ;;
        --structs_dir)
            structs_dir="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
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

if [ -z "$pae_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert pae_dir to absolute path if relative
if [[ "$pae_dir" != /* ]]; then
    pae_dir="$(cd "$pae_dir" && pwd)"
fi
# Convert structs_dir to absolute path if relative
if [[ -n "$structs_dir" && "$structs_dir" != /* ]]; then
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
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/structure_metrics

args=(--pae_dir "$pae_dir")
if [[ -n "$structs_dir" ]]; then
    args+=(--structs_dir "$structs_dir")
fi
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi

python run_global_confidence.py "${args[@]}"
