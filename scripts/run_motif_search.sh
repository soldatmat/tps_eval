#!/bin/bash

USAGE="--fasta_path <fasta_path> [<motif1> <motif2> ...]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path            Path to the FASTA file (required)"
    echo "  <motif1> <motif2> ...   Any number of motifs to search for in the sequences"
    echo "  -h, --help              Show this help message and exit"
    echo
}

# Parse long options manually
motifs=()
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
            motifs+=("$key")
            shift
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

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/../src/sequence_metrics"

julia run_motif_search.jl "$fasta_path" "${motifs[@]}"
