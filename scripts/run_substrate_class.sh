#!/bin/bash

USAGE="[--sequence_topk <csv>] [--embedding_topk <csv>] [--structural_topk <csv>] --label_file <csv> --calibration <json> [--pocket_csv <csv>] [--ee_csv <csv>] [--top_k <N>] --out <path>"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Substrate-class combiner. Fuses the label-agnostic k-NN substrate vote (the three"
    echo "similarity-space top-k CSVs + the substrate label file + the substrate calibration)"
    echo "with the pocket_descriptors catalytic_pocket_volume band and the EnzymeExplorer"
    echo "sequence-only per-substrate signal into <input>_substrate_class.csv."
    echo
    echo "Arguments:"
    echo "  --sequence_topk     <input>_max_sequence_identity_topk.csv (score = identity %)"
    echo "  --embedding_topk    <input>_min_embedding_distance_topk.csv (score = L2 distance)"
    echo "  --structural_topk   <structs_dir>_structural_identity_topk.csv (score = TM-score)"
    echo "  --label_file        SUBSTRATE reference_id,label CSV (substrate_labels.csv)"
    echo "  --calibration       Substrate calibration JSON (knn_calibration_substrate.json)"
    echo "  --pocket_csv        <structs_dir>_pocket_descriptors.csv (optional)"
    echo "  --ee_csv            <input>_enzyme_explorer_sequence_only.csv (optional)"
    echo "  --top_k             Cap neighbours per query (default: all present)"
    echo "  --out               Output predictions CSV (keyed by ID)"
    echo "  -h, --help          Show this help message and exit"
    echo
    echo "At least one of the three --*_topk CSVs is required."
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi
case "$1" in
    -h|--help) Help; exit 0 ;;
esac

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
        --sequence_topk|--embedding_topk|--structural_topk|--label_file|--calibration|--pocket_csv|--ee_csv|--out)
            passthru+=("$1" "$(abspath "$2")"); shift 2 ;;
        --top_k)
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

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting substrate_class..."
python run_substrate_class.py "${passthru[@]}"
rc=$?
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished substrate_class."
# Propagate python's exit code so a failed run FAILS the SLURM job (else the
# orchestrator's afterok dependents run on missing output -- a trailing echo would
# otherwise mask the failure with exit 0).
exit $rc
