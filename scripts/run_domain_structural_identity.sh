#!/bin/bash

USAGE="--structs_dir <structs_dir> [--known_domain_structures_root <dir>] [--save_path <csv>] [--n_jobs <n>] [--n_iters <n>] [--keep_detected_domains <dir>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir                   Directory of generated structures; EE detects domains in the .pdb files (ID = stem) (required)"
    echo "  --known_domain_structures_root  Directory of known-TPS reference DOMAIN structures (optional; default"
    echo "                                  \$ENZYME_EXPLORER_PATH/data/detected_domains/martsDB_detected_domains/domains)"
    echo "  --save_path                     Output CSV path (optional; default <structs_dir>_domain_structural_identity.csv)"
    echo "  --n_jobs                        Parallel jobs for detection (optional; default 10)"
    echo "  --n_iters                       EnzymeExplorer detection iterations (optional; default 3)"
    echo "  --keep_detected_domains         If given, keep the per-design detected domain .pdb files here (optional)"
    echo "  -h, --help                      Show this help message and exit"
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir) structs_dir="$2"; shift 2 ;;
        --known_domain_structures_root) known_domain_structures_root="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --n_jobs) n_jobs="$2"; shift 2 ;;
        --n_iters) n_iters="$2"; shift 2 ;;
        --keep_detected_domains) keep_detected_domains="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [ -z "$structs_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Absolute paths
[[ "$structs_dir" != /* ]] && structs_dir="$(cd "$structs_dir" && pwd)"
if [[ -n "$known_domain_structures_root" && "$known_domain_structures_root" != /* ]]; then
    known_domain_structures_root="$(cd "$known_domain_structures_root" && pwd)"
fi
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi
if [[ -n "$keep_detected_domains" && "$keep_detected_domains" != /* ]]; then
    mkdir -p "$keep_detected_domains"
    keep_detected_domains="$(cd "$keep_detected_domains" && pwd)"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load ENZYME_EXPLORER_ENV, ENZYME_EXPLORER_PATH

# Detection needs EnzymeExplorer's detect_domains AND foldseek (domain_alignment)
# in the SAME env. The EnzymeExplorer env (enzyme_explorer_prod on Aurum) has both.
eval "$(conda shell.bash hook)"
conda activate "$ENZYME_EXPLORER_ENV"
# Fix for compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29/.30
# (required by the EE env's pandas / PyMOL). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

# detect_domains imports EnzymeExplorer (enzymeexplorer.*) and the tool reads the
# template base PDBs from $ENZYME_EXPLORER_PATH/data/alphafold_structs. Export it.
export PYTHONPATH="$ENZYME_EXPLORER_PATH:${PYTHONPATH:-}"
export ENZYME_EXPLORER_PATH

# Default the reference DOMAINS root to EE's curated martsDB detected domains.
if [[ -z "$known_domain_structures_root" ]]; then
    known_domain_structures_root="$ENZYME_EXPLORER_PATH/data/detected_domains/martsDB_detected_domains/domains"
fi

cd src/structure_metrics

args=("$structs_dir" "$known_domain_structures_root")
[[ -n "$save_path" ]] && args+=(--save_path "$save_path")
[[ -n "$n_jobs" ]] && args+=(--n_jobs "$n_jobs")
[[ -n "$n_iters" ]] && args+=(--n_iters "$n_iters")
[[ -n "$keep_detected_domains" ]] && args+=(--keep_detected_domains "$keep_detected_domains")

python run_domain_structural_identity.py "${args[@]}"
