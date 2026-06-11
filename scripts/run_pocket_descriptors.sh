#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir   Directory of structures (.pdb/.cif) or AF3 af_output; file stem = ID (required)"
    echo "  --save_path     Output CSV path (optional; default <structs_dir>_pocket_descriptors.csv)"
    echo "  -h, --help      Show this help message and exit"
    echo
}

# Parse long options manually
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
. ./paths.sh # Load POCKET_ENV, P2RANK_PATH

eval "$(conda shell.bash hook)"
conda activate "$POCKET_ENV"
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
# P2Rank is the optional ML cross-check. Reference it via P2RANK_PATH (per-install
# absolute path in paths.sh, like SOLUPROT_PATH). The 'prank' launcher needs `java`
# on PATH, which the activated POCKET_ENV provides (openjdk). If P2RANK_PATH is unset
# the p2rank_* columns are NaN (fpocket still runs).
if [[ -n "$P2RANK_PATH" && -x "$P2RANK_PATH/prank" ]]; then
    args+=(--p2rank_bin "$P2RANK_PATH/prank")
else
    echo "[note] P2RANK_PATH unset or $P2RANK_PATH/prank not executable; skipping P2Rank cross-check."
fi

python run_pocket_descriptors.py "${args[@]}"
