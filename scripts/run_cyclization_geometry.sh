#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>] [--aromatic_cutoff <A>] [--farchain_bonds <n>] [--ion_resnames <R...>] [--min_substrate_carbons <n>] [--substrate_resname <R>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir            Directory of structures (.pdb/.cif) or AF3 af_output; file stem = ID (required)"
    echo "  --save_path              Output CSV path (optional; default <structs_dir>_cyclization_geometry.csv)"
    echo "  --aromatic_cutoff        Substrate-carbon -> aromatic centroid cation-pi cutoff (A) (optional; default 6.0)"
    echo "  --farchain_bonds         Bonds from C1 to count a carbon as distal for fold-back (optional; default 6)"
    echo "  --ion_resnames           Ion HETATM residue names to exclude from substrate detection (optional; default MG MN)"
    echo "  --min_substrate_carbons  Min carbons (with >=1 P) to count as a prenyl-PP substrate (optional; default 5)"
    echo "  --substrate_resname      Force a specific ligand residue name as the substrate (optional; default auto-detect)"
    echo "  -h, --help               Show this help message and exit"
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir) structs_dir="$2"; shift 2 ;;
        --save_path) save_path="$2"; shift 2 ;;
        --aromatic_cutoff) aromatic_cutoff="$2"; shift 2 ;;
        --farchain_bonds) farchain_bonds="$2"; shift 2 ;;
        --min_substrate_carbons) min_substrate_carbons="$2"; shift 2 ;;
        --substrate_resname) substrate_resname="$2"; shift 2 ;;
        --ion_resnames)
            shift
            ion_resnames=()
            while [[ $# -gt 0 && "$1" != --* ]]; do ion_resnames+=("$1"); shift; done
            ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [ -z "$structs_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

if [[ "$structs_dir" != /* ]]; then
    structs_dir="$(cd "$structs_dir" && pwd)"
fi
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Karolina compute-node libstdc++/GLIBCXX_3.4.29 fix: prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c \"import sys, json; print(json.load(sys.stdin)['active_prefix_name'])\")"
echo "Using python: $(which python)"

cd src/structure_metrics

args=("$structs_dir")
if [[ -n "$save_path" ]]; then args+=(--save_path "$save_path"); fi
if [[ -n "$aromatic_cutoff" ]]; then args+=(--aromatic_cutoff "$aromatic_cutoff"); fi
if [[ -n "$farchain_bonds" ]]; then args+=(--farchain_bonds "$farchain_bonds"); fi
if [[ -n "$min_substrate_carbons" ]]; then args+=(--min_substrate_carbons "$min_substrate_carbons"); fi
if [[ -n "$substrate_resname" ]]; then args+=(--substrate_resname "$substrate_resname"); fi
if [[ ${#ion_resnames[@]} -gt 0 ]]; then args+=(--ion_resnames "${ion_resnames[@]}"); fi

python run_cyclization_geometry.py "${args[@]}"
