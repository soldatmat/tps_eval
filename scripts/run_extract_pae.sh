#!/bin/bash

# Cluster-agnostic wrapper for src/alphafold/extract_pae.py: extract the per-structure
# PAE matrices (+ pTM/iPTM) from an AlphaFold3 af_output tree into the shared
# <ID>_pae.npz schema consumed by interdomain_pae / global_confidence. Used as the
# post-fold step of the orchestrator's `--fold alphafold3` path (after all AF3 jobs).

USAGE="--structs_dir <dir> --pae_dir <out_dir> | --af_output <dir> --pae_dir <out_dir>"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir   Structs dir whose sibling/child af_output/ holds the AF3 trees"
    echo "  --af_output     AF3 af_output dir directly (alternative to --structs_dir)"
    echo "  --pae_dir       Output dir for <ID>_pae.npz files (required)"
    echo "  -h, --help      Show this help message and exit"
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi

declare -a passthru=()
abspath() {
    local p="$1"
    if [[ "$p" != /* ]]; then p="$(cd "$(dirname "$p")" 2>/dev/null && pwd)/$(basename "$p")"; fi
    echo "$p"
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --structs_dir|--af_output|--pae_dir) passthru+=("$1" "$(abspath "$2")"); shift 2 ;;
        --no-skip_existing) passthru+=("$1"); shift ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/alphafold

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting AlphaFold3 PAE..."
python extract_pae.py "${passthru[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished extract_pae."
