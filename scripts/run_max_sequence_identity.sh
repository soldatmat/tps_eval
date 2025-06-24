#!/bin/bash

USAGE="--fasta_path <fasta_path> [--train_path <train_path>]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "If reference sequences are not provided, the script will use the sequences from the FASTA file itself"
    echo "and return the second maximum sequence identity \(first will be 100% with itself\)."
    echo
    echo "Arguments:"
    echo "  --fasta_path   Path to the FASTA file (required)"
    echo "  --train_path   Path to the reference FASTA file (optional)"
    echo "  -h, --help     Show this help message and exit"
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
        --train_path)
            train_path="$2"
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


if [[ -z "$fasta_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert fasta_path to absolute path if it's relative
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi

# Convert train_path to absolute path if it's set and relative
if [[ -n "$train_path" ]] && [[ "$train_path" != /* ]]; then
    train_path="$(cd "$(dirname "$train_path")" && pwd)/$(basename "$train_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/../src/sequence_metrics"

if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    julia run_max_sequence_identity.jl "$fasta_path" "$train_path"
else
    julia run_max_sequence_identity.jl "$fasta_path"
fi
