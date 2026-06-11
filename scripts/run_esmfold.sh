#!/bin/bash

USAGE="--fasta_path <fasta_path> --save_dir <save_dir> [--no-skip_existing] [--chunk_size <n>] [--device <cuda|cpu>] [--pae_dir <dir>] [--no-save_pae]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path        Path to the FASTA file of sequences to fold (required)"
    echo "  --save_dir          Directory to write <ID>.pdb structures into (required)"
    echo "  --no-skip_existing  Re-fold sequences even if <ID>.pdb already exists (optional)"
    echo "  --chunk_size        Force ESMFold trunk chunk size, lower = less GPU memory (optional)"
    echo "  --device            Torch device cuda/cpu (optional; default cuda if available)"
    echo "  --pae_dir           Directory for <ID>_pae.npz PAE matrices (optional; default <save_dir>_pae/)"
    echo "  --no-save_pae       Do not save PAE matrices (optional; default: save them)"
    echo "  -h, --help          Show this help message and exit"
    echo
}

# Parse long options manually
skip_existing_flag=""
save_pae_flag=""
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path)
            fasta_path="$2"
            shift 2
            ;;
        --save_dir)
            save_dir="$2"
            shift 2
            ;;
        --no-skip_existing)
            skip_existing_flag="--no-skip_existing"
            shift
            ;;
        --chunk_size)
            chunk_size="$2"
            shift 2
            ;;
        --device)
            device="$2"
            shift 2
            ;;
        --pae_dir)
            pae_dir="$2"
            shift 2
            ;;
        --no-save_pae)
            save_pae_flag="--no-save_pae"
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

if [ -z "$fasta_path" ] || [ -z "$save_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert fasta_path to absolute path if relative
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi
# Convert save_dir to absolute path if relative (create it first so cd/dirname works)
mkdir -p "$save_dir"
if [[ "$save_dir" != /* ]]; then
    save_dir="$(cd "$save_dir" && pwd)"
fi
# Convert pae_dir to absolute path if relative (create it first so cd/dirname works)
if [[ -n "$pae_dir" ]]; then
    mkdir -p "$pae_dir"
    if [[ "$pae_dir" != /* ]]; then
        pae_dir="$(cd "$pae_dir" && pwd)"
    fi
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load ESMFOLD_ENV

eval "$(conda shell.bash hook)"
conda activate "$ESMFOLD_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/esmfold

args=("$fasta_path" --save_dir "$save_dir")
if [[ -n "$skip_existing_flag" ]]; then
    args+=("$skip_existing_flag")
fi
if [[ -n "$chunk_size" ]]; then
    args+=(--chunk_size "$chunk_size")
fi
if [[ -n "$device" ]]; then
    args+=(--device "$device")
fi
if [[ -n "$pae_dir" ]]; then
    args+=(--pae_dir "$pae_dir")
fi
if [[ -n "$save_pae_flag" ]]; then
    args+=("$save_pae_flag")
fi

python run_esmfold.py "${args[@]}"
