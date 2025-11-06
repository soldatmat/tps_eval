#!/bin/bash

USAGE="--fasta_path <fasta_path>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path    Path to the FASTA file (required)"
    echo "  -h, --help      Show this help message and exit"
    echo
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path)
            fasta_path="$2"
            shift
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

python src/esm/extract_embeddings.py \
    esm1b_t33_650M_UR50S \
    "$fasta_path" \
    --repr_layers 33 \
    --include mean
