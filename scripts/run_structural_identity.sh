#!/bin/bash

USAGE="--structs_dir <structs_dir> --known_structs_dir <known_structs_dir> [--save_path <save_path> --top_k <N>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir         Directory of generated structures (.pdb/.cif) (required)"
    echo "  --known_structs_dir   Directory of known-TPS reference structures (required)"
    echo "  --save_path           Output CSV path (optional; default <structs_dir>_structural_identity.csv)"
    echo "  --top_k               If >=1, also write <structs_dir>_structural_identity_topk.csv (query_id,rank,neighbour_id,score; score = TM-score, LARGER closer)"
    echo "  -h, --help            Show this help message and exit"
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir) structs_dir="$2"; shift 2 ;;
        --known_structs_dir) known_structs_dir="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --top_k) top_k="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$structs_dir" || -z "$known_structs_dir" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Absolute paths
[[ "$structs_dir" != /* ]] && structs_dir="$(cd "$structs_dir" && pwd)"
[[ "$known_structs_dir" != /* ]] && known_structs_dir="$(cd "$known_structs_dir" && pwd)"
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV (provides foldseek)

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/structure_metrics

args=("$structs_dir" "$known_structs_dir")
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi
if [[ -n "$top_k" ]]; then
    args+=(--top_k "$top_k")
fi

python run_structural_identity.py "${args[@]}"
