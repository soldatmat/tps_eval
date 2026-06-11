#!/bin/bash

USAGE="<calibrate|predict> [--sequence_topk <csv>] [--embedding_topk <csv>] [--structural_topk <csv>] --label_file <csv> --out <path> [--calibration <json>] [--labeling <name>] [--top_k <N>] [--target_accuracy <a>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Label-agnostic k-NN coarse-label transfer. Distance-weighted vote of nearest"
    echo "MARTS-DB neighbours per design, ensembled across the three tps_eval similarity"
    echo "spaces, with leave-one-out calibration on MARTS-DB."
    echo
    echo "Subcommands:"
    echo "  calibrate   LOO calibration on MARTS-DB SELF top-k CSVs -> --out JSON artifact."
    echo "  predict     Transfer labels to designs from their top-k CSVs + --calibration JSON."
    echo
    echo "Arguments:"
    echo "  --sequence_topk     <input>_max_sequence_identity_topk.csv (score = identity %)"
    echo "  --embedding_topk    <input>_min_embedding_distance_topk.csv (score = L2 distance)"
    echo "  --structural_topk   <structs_dir>_structural_identity_topk.csv (score = TM-score)"
    echo "  --label_file        CSV mapping reference_id,label (the labeling is the INPUT)"
    echo "  --out               Output path (calibration JSON, or predictions CSV)"
    echo "  --calibration       Calibration JSON (predict mode only)"
    echo "  --labeling          Name recorded in the artifact (calibrate mode; e.g. first_cyclization)"
    echo "  --top_k             Cap neighbours per query (default: all present)"
    echo "  --target_accuracy   Accuracy floor for tau selection (calibrate mode; default 0.5)"
    echo "  -h, --help          Show this help message and exit"
    echo
    echo "At least one of the three --*_topk CSVs is required."
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi
case "$1" in
    -h|--help) Help; exit 0 ;;
esac

cmd="$1"; shift
if [[ "$cmd" != "calibrate" && "$cmd" != "predict" ]]; then
    echo "First argument must be 'calibrate' or 'predict'."; Help; exit 1
fi

# Pass remaining args through; convert relative path-bearing args to absolute.
declare -a passthru=()
abspath() {
    local p="$1"
    if [[ "$p" != /* ]]; then
        if [[ -e "$p" ]]; then p="$(cd "$(dirname "$p")" && pwd)/$(basename "$p")";
        else p="$(pwd)/$p"; fi
    fi
    echo "$p"
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sequence_topk|--embedding_topk|--structural_topk|--label_file|--out|--calibration)
            passthru+=("$1" "$(abspath "$2")"); shift 2 ;;
        --labeling|--top_k|--target_accuracy)
            passthru+=("$1" "$2"); shift 2 ;;
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
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/knn

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting knn_label_transfer ($cmd)..."
python run_knn_label_transfer.py "$cmd" "${passthru[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished knn_label_transfer ($cmd)."
