#!/bin/bash

USAGE="--fasta_path <fasta_path> --num_seqs <num_seqs> --save_to <save_to> [--return_counts True|False]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path                Path to the FASTA file (required)"
    echo "  --num_seqs                  Number of sequences to sample (required)"
    echo "  --save_to                   Path to save the sampled lengths (required)"
    echo "  --return_counts             Whether to return unique length with their counts (default: False)"
    echo "  -h, --help                  Show this help message and exit"
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
        --num_seqs)
            num_seqs="$2"
            shift
            shift
            ;;
        --save_to)
            save_to="$2"
            shift
            shift
            ;;
        --return_counts)
            return_counts="$2"
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


if [[ -z "$fasta_path" ]] || [[ -z "$num_seqs" ]] || [[ -z "$save_to" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

if [[ -z "$return_counts" ]]; then
    return_counts="False"  # Default value
fi

# Convert fasta_path to absolute path if it's relative
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/.."

eval "$(conda shell.bash hook)"
conda activate terpene_generation

python "$SCRIPT_DIR/../src/sample_length/sample_length.py" \
    --fasta_path "$fasta_path" \
    --num_seqs "$num_seqs" \
    --save_to "$save_to" \
    --return_counts "$return_counts"
