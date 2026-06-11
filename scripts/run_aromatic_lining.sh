#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>] [--cutoff <A>] [--cation_pi_min <A>] [--cation_pi_max <A>] [--face_angle_deg <deg>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir     Directory of structures (.pdb/.cif); file stem = ID (required)"
    echo "  --save_path       Output CSV path (optional; default <structs_dir>_aromatic_lining.csv)"
    echo "  --cutoff          Pocket shell radius in A around the metal point (optional; default 10)"
    echo "  --cation_pi_min   Min ring-centroid->locus distance for inward hit (optional; default 3.5)"
    echo "  --cation_pi_max   Max ring-centroid->locus distance for inward hit (optional; default 6.0)"
    echo "  --face_angle_deg  Max ring-normal vs locus angle for face-on hit (optional; default 45)"
    echo "  -h, --help        Show this help message and exit"
    echo
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir)
            structs_dir="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
            shift 2
            ;;
        --cutoff)
            cutoff="$2"
            shift 2
            ;;
        --cation_pi_min)
            cation_pi_min="$2"
            shift 2
            ;;
        --cation_pi_max)
            cation_pi_max="$2"
            shift 2
            ;;
        --face_angle_deg)
            face_angle_deg="$2"
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

if [ -z "$structs_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert structs_dir to absolute path if relative
if [[ "$structs_dir" != /* ]]; then
    structs_dir="$(cd "$structs_dir" && pwd)"
fi
# Convert save_path to absolute path if relative
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

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

cd src/structure_metrics

args=("$structs_dir")
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi
if [[ -n "$cutoff" ]]; then
    args+=(--cutoff "$cutoff")
fi
if [[ -n "$cation_pi_min" ]]; then
    args+=(--cation_pi_min "$cation_pi_min")
fi
if [[ -n "$cation_pi_max" ]]; then
    args+=(--cation_pi_max "$cation_pi_max")
fi
if [[ -n "$face_angle_deg" ]]; then
    args+=(--face_angle_deg "$face_angle_deg")
fi

python run_aromatic_lining.py "${args[@]}"
