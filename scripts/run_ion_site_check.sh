#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>] [--site_radius <A>] [--coord_cutoff <A>] [--min_coord_contacts <n>] [--ion_resnames <R...>] [--diphosphate_resnames <R...>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir          Directory of structures (.pdb/.cif); file stem = ID (required)"
    echo "  --save_path            Output CSV path (optional; default <structs_dir>_ion_site_check.csv)"
    echo "  --site_radius          In-site distance (A) from the cage centroid (optional; default 5.0)"
    echo "  --coord_cutoff         Mg-O coordination cutoff (A) (optional; default 2.8)"
    echo "  --min_coord_contacts   Min coordinating-O contacts for well-placed (optional; default 2)"
    echo "  --ion_resnames         Ion HETATM residue names (optional; default MG MN)"
    echo "  --diphosphate_resnames Diphosphate HETATM residue names (optional; default POP PPV PPK)"
    echo "  -h, --help             Show this help message and exit"
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
        --site_radius)
            site_radius="$2"
            shift 2
            ;;
        --coord_cutoff)
            coord_cutoff="$2"
            shift 2
            ;;
        --min_coord_contacts)
            min_coord_contacts="$2"
            shift 2
            ;;
        --ion_resnames)
            shift
            ion_resnames=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                ion_resnames+=("$1")
                shift
            done
            ;;
        --diphosphate_resnames)
            shift
            diphosphate_resnames=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                diphosphate_resnames+=("$1")
                shift
            done
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
if [[ -n "$site_radius" ]]; then
    args+=(--site_radius "$site_radius")
fi
if [[ -n "$coord_cutoff" ]]; then
    args+=(--coord_cutoff "$coord_cutoff")
fi
if [[ -n "$min_coord_contacts" ]]; then
    args+=(--min_coord_contacts "$min_coord_contacts")
fi
if [[ ${#ion_resnames[@]} -gt 0 ]]; then
    args+=(--ion_resnames "${ion_resnames[@]}")
fi
if [[ ${#diphosphate_resnames[@]} -gt 0 ]]; then
    args+=(--diphosphate_resnames "${diphosphate_resnames[@]}")
fi

python run_ion_site_check.py "${args[@]}"
