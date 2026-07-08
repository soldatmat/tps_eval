#!/bin/bash

USAGE="--fasta_path <fasta_path> [--target_substrate <CODE>] [--smiles <SMILES>] [--device cuda|cpu] [--model_dpath <dir>] [--batch_size <n>] [--out_suffix <s>] [--save_path <path>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path        Path to the FASTA file (required)"
    echo "  --target_substrate  Substrate code (e.g. GPP, FPP, GGPP, GFPP) resolved to a SMILES"
    echo "  --smiles            Explicit substrate SMILES (overrides --target_substrate)"
    echo "  --device            Torch device (default cuda; use cpu on CPU nodes)"
    echo "  --model_dpath       CataPro model dir (default vendor/CataPro/models)"
    echo "  --batch_size        CataPro batch size (default 64)"
    echo "  --out_suffix        Output suffix -> <fasta>_<suffix>.csv (default catapro)"
    echo "  --save_path         Output CSV path (optional)"
    echo "  -h, --help          Show this help message and exit"
    echo
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path)        fasta_path="$2"; shift 2 ;;
        --target_substrate)  target_substrate="$2"; shift 2 ;;
        --smiles)            smiles="$2"; shift 2 ;;
        --device)            device="$2"; shift 2 ;;
        --model_dpath)       model_dpath="$2"; shift 2 ;;
        --batch_size)        batch_size="$2"; shift 2 ;;
        --out_suffix)        out_suffix="$2"; shift 2 ;;
        --save_path)         save_path="$2"; shift 2 ;;
        -h|--help)           Help; exit 0 ;;
        *)                   echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [ -z "$fasta_path" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert fasta_path to absolute path if it's relative
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
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
. ./paths.sh # Load CATAPRO_ENV

eval "$(conda shell.bash hook)"
conda activate "$CATAPRO_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/sequence_metrics

args=("$fasta_path")
if [[ -n "$target_substrate" ]]; then args+=(--target_substrate "$target_substrate"); fi
if [[ -n "$smiles" ]]; then args+=(--smiles "$smiles"); fi
if [[ -n "$device" ]]; then args+=(--device "$device"); fi
if [[ -n "$model_dpath" ]]; then args+=(--model_dpath "$model_dpath"); fi
if [[ -n "$batch_size" ]]; then args+=(--batch_size "$batch_size"); fi
if [[ -n "$out_suffix" ]]; then args+=(--out_suffix "$out_suffix"); fi
if [[ -n "$save_path" ]]; then args+=(--save_path "$save_path"); fi

python run_catapro.py "${args[@]}"
