#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <csv> --top_n <n> --max_seqs <n>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Broad STRUCTURE homology search of designs vs AlphaFold-Swiss-Prot (foldseek)."
    echo "Reports each design's top hit across all proteins and whether the top-N hits"
    echo "are terpene synthases. Output: <structs_dir>_foldseek_swissprot_search.csv."
    echo
    echo "Arguments:"
    echo "  --structs_dir  Directory of generated structures (af_output or flat dir) (required)"
    echo "  --save_path    Output CSV path (optional; default <structs_dir>_foldseek_swissprot_search.csv)"
    echo "  --top_n        Top-N hits per query (optional; default 25)"
    echo "  --max_seqs     foldseek --max-seqs prefilter (optional; default 300)"
    echo "  -h, --help     Show this help message and exit"
    echo
    echo "Requires paths.sh: AFDB_SWISSPROT_DB (foldseek AFDB-Swiss-Prot DB) and TPS_ACCESSIONS."
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir) structs_dir="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --top_n) top_n="$2"; shift 2 ;;
        --max_seqs) max_seqs="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$structs_dir" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Absolute paths
[[ "$structs_dir" != /* ]] && structs_dir="$(cd "$structs_dir" && pwd)"
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV (provides foldseek), AFDB_SWISSPROT_DB, TPS_ACCESSIONS

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

if [[ -z "$AFDB_SWISSPROT_DB" ]]; then
    echo "Error: AFDB_SWISSPROT_DB is not set in paths.sh."
    exit 1
fi
TPS_ACCESSIONS="${TPS_ACCESSIONS:-$(pwd)/src/homology_search/tps_uniprot_accessions.txt}"

cd src/homology_search

args=("$structs_dir" "$AFDB_SWISSPROT_DB" "$TPS_ACCESSIONS")
[[ -n "$save_path" ]] && args+=(--save_path "$save_path")
[[ -n "$top_n" ]] && args+=(--top_n "$top_n")
[[ -n "$max_seqs" ]] && args+=(--max_seqs "$max_seqs")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting foldseek_swissprot_search..."
python run_foldseek_swissprot_search.py "${args[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished foldseek_swissprot_search."
