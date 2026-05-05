#!/bin/bash

USAGE="--structures_selection_csv <csv> --structures_root <dir> --known_structures_root <dir> --output_root <dir>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "For each row in the selection CSV, render PyMOL images that color the new structure"
    echo "by sequence-similarity (BLOSUM90) to the matched known structure, plus a side-by-side alignment."
    echo
    echo "Arguments:"
    echo "  --structures_selection_csv          CSV with one row per (query, target) pair (required)"
    echo "  --structures_column_name            Column with new-structure name (optional, default 'query')"
    echo "  --known_structures_column_name      Column with matched known-structure name (optional, default 'max_alntmscore_target')"
    echo "  --similarity_metric_name            Metric name used for output dir naming (optional, default 'similarity')"
    echo "  --structures_root                   Directory holding new structure PDBs (required)"
    echo "  --known_structures_root             Directory holding known reference PDBs (required)"
    echo "  --output_root                       Directory to write PNGs and PyMOL sessions (required)"
    echo "  --no-store_pymol_sessions           Skip saving .pse session files (optional)"
    echo "  --no-rerun_existing                 Skip jobs whose output already exists (optional)"
    echo "  -h, --help                          Show this help message and exit"
    echo
}

# Collect extra arguments for plot_residue_similarity.py
extra_args=()

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structures_selection_csv)
            structures_selection_csv="$2"
            shift 2
            ;;
        --structures_column_name)
            extra_args+=(--structures_column_name "$2")
            shift 2
            ;;
        --known_structures_column_name)
            extra_args+=(--known_structures_column_name "$2")
            shift 2
            ;;
        --similarity_metric_name)
            extra_args+=(--similarity_metric_name "$2")
            shift 2
            ;;
        --structures_root)
            structures_root="$2"
            shift 2
            ;;
        --known_structures_root)
            known_structures_root="$2"
            shift 2
            ;;
        --output_root)
            output_root="$2"
            shift 2
            ;;
        --store_pymol_sessions)
            extra_args+=(--store_pymol_sessions)
            shift
            ;;
        --no-store_pymol_sessions)
            extra_args+=(--no-store_pymol_sessions)
            shift
            ;;
        --rerun_existing)
            extra_args+=(--rerun_existing)
            shift
            ;;
        --no-rerun_existing)
            extra_args+=(--no-rerun_existing)
            shift
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

if [[ -z "$structures_selection_csv" || -z "$structures_root" || -z "$known_structures_root" || -z "$output_root" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert paths to absolute if relative
for var in structures_selection_csv structures_root known_structures_root output_root; do
    val="${!var}"
    if [[ "$val" != /* ]]; then
        if [[ -e "$val" ]]; then
            printf -v "$var" "%s" "$(cd "$(dirname "$val")" && pwd)/$(basename "$val")"
        else
            printf -v "$var" "%s" "$(pwd)/$val"
        fi
    fi
done

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting plot_residue_similarity..."
python -m src.pymol.plot_residue_similarity \
    --structures_selection_csv "$structures_selection_csv" \
    --structures_root "$structures_root" \
    --known_structures_root "$known_structures_root" \
    --output_root "$output_root" \
    "${extra_args[@]}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished plot_residue_similarity."
