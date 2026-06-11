#!/bin/bash

USAGE="--structs_dir <structs_dir> --pae_dir <pae_dir> [--save_path <save_path>] [--detections_json <json>] [--per_pair] [--n_jobs <n>] [--n_iters <n>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir       Directory of generated structures; EE detects domains in the .pdb files (ID = stem) (required)"
    echo "  --pae_dir           Directory of <ID>_pae.npz PAE matrices saved at fold time (required)"
    echo "  --save_path         Output CSV path (optional; default <structs_dir>_interdomain_pae.csv)"
    echo "  --detections_json   Existing EE domain-detection JSON sidecar to PARSE instead of re-detecting (optional)"
    echo "  --per_pair          Also emit a pae_<A>_<B> column per inter-domain pair (optional)"
    echo "  --n_jobs            Parallel jobs for detection (optional; default 10)"
    echo "  --n_iters           EnzymeExplorer detection iterations (optional; default 3)"
    echo "  -h, --help          Show this help message and exit"
    echo
}

# Parse long options manually
per_pair_flag=""
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir)
            structs_dir="$2"
            shift 2
            ;;
        --pae_dir)
            pae_dir="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
            shift 2
            ;;
        --detections_json)
            detections_json="$2"
            shift 2
            ;;
        --per_pair)
            per_pair_flag="--per_pair"
            shift
            ;;
        --n_jobs)
            n_jobs="$2"
            shift 2
            ;;
        --n_iters)
            n_iters="$2"
            shift 2
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

if [ -z "$structs_dir" ] || [ -z "$pae_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert structs_dir to absolute path if relative
if [[ "$structs_dir" != /* ]]; then
    structs_dir="$(cd "$structs_dir" && pwd)"
fi
# Convert pae_dir to absolute path if relative
if [[ "$pae_dir" != /* ]]; then
    pae_dir="$(cd "$pae_dir" && pwd)"
fi
# Convert save_path to absolute path if relative
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi
# Convert detections_json to absolute path if relative (parent must exist)
if [[ -n "$detections_json" && "$detections_json" != /* ]]; then
    detections_json="$(cd "$(dirname "$detections_json")" && pwd)/$(basename "$detections_json")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load ENZYME_EXPLORER_ENV, ENZYME_EXPLORER_PATH variables

eval "$(conda shell.bash hook)"
conda activate "$ENZYME_EXPLORER_ENV"
# Fix for compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29/.30
# (required by the EE env's pandas / PyMOL). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

# interdomain_pae imports EnzymeExplorer (enzymeexplorer.*) via domain_composition,
# so run with the EE repo on PYTHONPATH (mirrors run_domain_composition.sh).
export PYTHONPATH="$ENZYME_EXPLORER_PATH:${PYTHONPATH:-}"

# Run the python entry from the tool dir (so `import interdomain_pae` resolves).
cd src/structure_metrics

args=("$structs_dir" --pae_dir "$pae_dir")
[[ -n "$save_path" ]] && args+=(--save_path "$save_path")
[[ -n "$detections_json" ]] && args+=(--detections_json "$detections_json")
[[ -n "$per_pair_flag" ]] && args+=("$per_pair_flag")
[[ -n "$n_jobs" ]] && args+=(--n_jobs "$n_jobs")
[[ -n "$n_iters" ]] && args+=(--n_iters "$n_iters")

python run_interdomain_pae.py "${args[@]}"
