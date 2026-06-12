#!/bin/bash

USAGE="--structures_selection_csv <csv> --structures_root <dir> --domain_structures_root <dir> --domains_pkl <pkl> --output_root <dir>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Render per-design PyMOL images of detected domains overlaid on the full structure."
    echo
    echo "Arguments:"
    echo "  --structures_selection_csv      CSV with one row per structure (required)"
    echo "  --structures_column_name        Column with structure name (optional, default 'ID')"
    echo "  --structures_file_suffix        Suffix appended to structure name when looking up <structures_root>/<name><suffix>.pdb (optional)"
    echo "  --structures_root               Directory holding full structure PDBs (required)"
    echo "  --domain_structures_root        Directory holding per-domain PDBs (required)"
    echo "  --domains_pkl                   Pickle from EnzymeExplorer with detected domains (required)"
    echo "  --output_root                   Directory to write PNGs and PyMOL sessions (required)"
    echo "  --no-store_pymol_sessions       Skip saving .pse session files (optional)"
    echo "  --no-rerun_existing             Skip jobs whose output already exists (optional)"
    echo "  -h, --help                      Show this help message and exit"
    echo
}

# Collect extra arguments for plot_domains.py
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
        --structures_file_suffix)
            extra_args+=(--structures_file_suffix "$2")
            shift 2
            ;;
        --structures_root)
            structures_root="$2"
            shift 2
            ;;
        --domain_structures_root)
            domain_structures_root="$2"
            shift 2
            ;;
        --domains_pkl)
            domains_pkl="$2"
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

if [[ -z "$structures_selection_csv" || -z "$structures_root" || -z "$domain_structures_root" || -z "$domains_pkl" || -z "$output_root" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert paths to absolute if relative
for var in structures_selection_csv structures_root domain_structures_root domains_pkl output_root; do
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
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting plot_domains..."
python -m src.pymol.plot_domains \
    --structures_selection_csv "$structures_selection_csv" \
    --structures_root "$structures_root" \
    --domain_structures_root "$domain_structures_root" \
    --domains_pkl "$domains_pkl" \
    --output_root "$output_root" \
    "${extra_args[@]}"
rc=$?
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished plot_domains."
# Propagate python's exit code so a failed render FAILS the SLURM job (a trailing
# echo would otherwise mask the failure with exit 0).
exit $rc
