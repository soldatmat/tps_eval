#!/bin/bash

USAGE="--fasta_path <fasta_path> [--train_path <ref.fasta> --backend <mmseqs2|diamond> --top_k <N> --threads <n> --sensitivity <s> --save_path <csv> --topk_save_path <csv>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Fast LOCAL (BLAST-style) sequence identity/similarity search, backend-pluggable"
    echo "over MMseqs2 (default) and DIAMOND. The LOCAL counterpart of max_sequence_identity"
    echo "(the GLOBAL full-length novelty metric, which is unaffected)."
    echo
    echo "If --train_path is omitted, self mode is used (each query's best hit / neighbours"
    echo "exclude itself)."
    echo
    echo "Arguments:"
    echo "  --fasta_path      Path to the query FASTA file (required)"
    echo "  --train_path      Path to the reference FASTA file (optional; omit for self mode)"
    echo "  --backend         mmseqs2 (default) or diamond"
    echo "  --top_k           If >=1, also write <input>_local_sequence_search_topk.csv"
    echo "                    (query_id,rank,neighbour_id,score; score = identity percent, LARGER closer)"
    echo "  --threads         Backend threads (optional; default \$SLURM_CPUS_PER_TASK or 4)"
    echo "  --sensitivity     Backend sensitivity knob (mmseqs2 -s value; diamond flag name)"
    echo "  --save_path       Metric CSV path (optional)"
    echo "  --topk_save_path  Top-k CSV path (optional)"
    echo "  -h, --help        Show this help message and exit"
    echo
    echo "Output: <input>_local_sequence_search.csv keyed by ID. Per-backend mapping:"
    echo "  local_sequence_identity   = mmseqs2 fident*100 | diamond pident"
    echo "  local_sequence_similarity = mmseqs2 NaN (no positives field) | diamond ppos"
    echo "  local_coverage            = mmseqs2 qcov*100 | diamond qcovhsp"
    echo
    echo "Requires paths.sh: TPS_EVAL_ENV (must contain mmseqs2 + diamond)."
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path) fasta_path="$2"; shift 2 ;;
        --train_path) train_path="$2"; shift 2 ;;
        --backend) backend="$2"; shift 2 ;;
        --top_k) top_k="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --sensitivity) sensitivity="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --topk_save_path) topk_save_path="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$fasta_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert paths to absolute when relative.
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi
if [[ -n "$train_path" && "$train_path" != /* ]]; then
    train_path="$(cd "$(dirname "$train_path")" && pwd)/$(basename "$train_path")"
fi
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi
if [[ -n "$topk_save_path" && "$topk_save_path" != /* ]]; then
    topk_save_path="$(cd "$(dirname "$topk_save_path")" && pwd)/$(basename "$topk_save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV (or LOCAL_SEQUENCE_SEARCH_ENV if a dedicated env was made)

# Allow a dedicated env override (set LOCAL_SEQUENCE_SEARCH_ENV in paths.sh if mmseqs2
# could not be co-installed into TPS_EVAL_ENV). Defaults to TPS_EVAL_ENV.
SEARCH_ENV="${LOCAL_SEQUENCE_SEARCH_ENV:-$TPS_EVAL_ENV}"

eval "$(conda shell.bash hook)"
conda activate "$SEARCH_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

# Default backend threads to the SLURM allocation when available.
threads="${threads:-${SLURM_CPUS_PER_TASK:-4}}"

cd src/sequence_metrics

args=("$fasta_path")
[[ -n "$train_path" ]] && args+=("$train_path")
args+=(--threads "$threads")
[[ -n "$backend" ]] && args+=(--backend "$backend")
[[ -n "$top_k" ]] && args+=(--top_k "$top_k")
[[ -n "$sensitivity" ]] && args+=(--sensitivity "$sensitivity")
[[ -n "$save_path" ]] && args+=(--save_path "$save_path")
[[ -n "$topk_save_path" ]] && args+=(--topk_save_path "$topk_save_path")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting local_sequence_search (backend=${backend:-mmseqs2})..."
python run_local_sequence_search.py "${args[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished local_sequence_search."
