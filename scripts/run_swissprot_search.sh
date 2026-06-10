#!/bin/bash

USAGE="--fasta_path <fasta_path> [--save_path <csv> --top_n <n> --threads <n> --sensitivity <flag>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Broad SEQUENCE homology search of designs vs Swiss-Prot (DIAMOND blastp)."
    echo "Reports each design's top hit across all proteins and whether the top-N hits"
    echo "are terpene synthases. Output: <fasta>_swissprot_search.csv (keyed by ID)."
    echo
    echo "Arguments:"
    echo "  --fasta_path   Path to the design FASTA file (required)"
    echo "  --save_path    Output CSV path (optional; default <fasta>_swissprot_search.csv)"
    echo "  --top_n        Top-N hits per query (optional; default 25)"
    echo "  --threads      DIAMOND threads (optional; default \$SLURM_CPUS_PER_TASK or 4)"
    echo "  --sensitivity  DIAMOND sensitivity flag (optional; default very-sensitive)"
    echo "  -h, --help     Show this help message and exit"
    echo
    echo "Requires paths.sh: SWISSPROT_DIAMOND_DB (DIAMOND DB) and TPS_ACCESSIONS."
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path) fasta_path="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --top_n) top_n="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --sensitivity) sensitivity="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$fasta_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Absolute paths
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV, SWISSPROT_DIAMOND_DB, TPS_ACCESSIONS

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

if [[ -z "$SWISSPROT_DIAMOND_DB" ]]; then
    echo "Error: SWISSPROT_DIAMOND_DB is not set in paths.sh."
    exit 1
fi
# TPS_ACCESSIONS defaults to the committed reference file if unset.
TPS_ACCESSIONS="${TPS_ACCESSIONS:-$(pwd)/src/homology_search/tps_uniprot_accessions.txt}"

# Default DIAMOND threads to the SLURM allocation when available.
threads="${threads:-${SLURM_CPUS_PER_TASK:-4}}"

cd src/homology_search

args=("$fasta_path" "$SWISSPROT_DIAMOND_DB" "$TPS_ACCESSIONS" --threads "$threads")
[[ -n "$save_path" ]] && args+=(--save_path "$save_path")
[[ -n "$top_n" ]] && args+=(--top_n "$top_n")
[[ -n "$sensitivity" ]] && args+=(--sensitivity "$sensitivity")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting swissprot_search (DIAMOND blastp)..."
python run_swissprot_search.py "${args[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished swissprot_search."
