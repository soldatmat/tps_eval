#!/bin/bash

USAGE="--fasta_path <fasta_path> [--save_path <save_path>] [--method swoop|masked] [--nogpu]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path    Path to the FASTA file (required)"
    echo "  --save_path     Output CSV path (optional; default <fasta>_esm_pseudo_perplexity.csv)"
    echo "  --method        swoop (fast single-pass approx, default) or masked (exact, slow)"
    echo "  --model_location ESM model name/path (optional; default esm1b_t33_650M_UR50S)"
    echo "  --nogpu         Do not use GPU even if available (optional)"
    echo "  -h, --help      Show this help message and exit"
    echo
}

nogpu_flag=""
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path)
            fasta_path="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
            shift 2
            ;;
        --method)
            method="$2"
            shift 2
            ;;
        --model_location)
            model_location="$2"
            shift 2
            ;;
        --nogpu)
            nogpu_flag="--nogpu"
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
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/sequence_metrics

args=("$fasta_path")
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi
if [[ -n "$method" ]]; then
    args+=(--method "$method")
fi
if [[ -n "$model_location" ]]; then
    args+=(--model_location "$model_location")
fi
if [[ -n "$nogpu_flag" ]]; then
    args+=("$nogpu_flag")
fi

python run_esm_pseudo_perplexity.py "${args[@]}"
