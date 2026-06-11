#!/bin/bash

USAGE="--structs_dir <structs_dir> --known_structs_dir <known_structs_dir> [--structural_topk <csv>] [--sequence_topk <csv>] [--sdr_panel <csv>] [--panel_cutoff <A>] [--map_tolerance <A>] [--tau_high <s>] [--tau_low <s>] [--save_path <path>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Specificity-determining-residue (SDR) divergence: flag each design that is"
    echo "GLOBALLY close to a known-product TPS but DIVERGES at the active-site residues"
    echo "that determine product specificity (the TEAS/HPS single-switch failure mode)."
    echo "Reuses the committed --top_k neighbour CSVs + Biopython structural superposition."
    echo
    echo "Arguments:"
    echo "  --structs_dir        Directory of design structures (.pdb/.cif or AF3 af_output) (required)"
    echo "  --known_structs_dir  Directory of known-TPS reference structures (required)"
    echo "  --structural_topk    <structs_dir>_structural_identity_topk.csv (TM-score; preferred)"
    echo "  --sequence_topk      <input>_max_sequence_identity_topk.csv (identity %; fallback)"
    echo "  --sdr_panel          Optional explicit SDR panel CSV (default: structure-derived)"
    echo "  --panel_cutoff       Structure-derived panel cutoff in A around the metal point (default 10)"
    echo "  --map_tolerance      Max Ca-Ca map distance in A after superposition (default 4)"
    echo "  --tau_high           Global-similarity floor in [0,1] (default 0.6)"
    echo "  --tau_low            SDR-identity ceiling in [0,1] (default 0.7)"
    echo "  --save_path          Output CSV path (default <structs_dir>_sdr_divergence.csv)"
    echo "  -h, --help           Show this help message and exit"
    echo
    echo "At least one of --structural_topk / --sequence_topk is required."
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi

abspath() {
    local p="$1"
    if [[ "$p" != /* ]]; then
        if [[ -e "$p" ]]; then p="$(cd "$(dirname "$p")" && pwd)/$(basename "$p")";
        else p="$(pwd)/$p"; fi
    fi
    echo "$p"
}

structs_dir=""
known_structs_dir=""
declare -a opts=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --structs_dir) structs_dir="$(abspath "$2")"; shift 2 ;;
        --known_structs_dir) known_structs_dir="$(abspath "$2")"; shift 2 ;;
        --structural_topk|--sequence_topk|--sdr_panel|--save_path)
            opts+=("$1" "$(abspath "$2")"); shift 2 ;;
        --panel_cutoff|--map_tolerance|--tau_high|--tau_low)
            opts+=("$1" "$2"); shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$structs_dir" || -z "$known_structs_dir" ]]; then
    echo "Usage: $0 $USAGE"; exit 1
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/specificity

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting sdr_divergence..."
python run_sdr_divergence.py "$structs_dir" "$known_structs_dir" "${opts[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished sdr_divergence."
