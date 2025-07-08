#!/bin/bash

USAGE="--embeddings_path <embeddings_path> [--train_embeddings_path <train_embeddings_path>]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "If reference sequences are not provided, the script will use the sequences from the FASTA file itself"
    echo "and return the second minimum embedding distance \(first will be 0. with itself\)."
    echo
    echo "Arguments:"
    echo "  --embeddings_path           Path to the CSV file with embeddings (required)"
    echo "  --train_embeddings_path     Path to the reference CSV file with embeddings (optional)"
    echo "  --train                     Turns on train data mode. "_self" results will be also copied as non-"_self" results."
    echo "  -h, --help                  Show this help message and exit"
    echo
}

# Parse long options manually
train_mode=false
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --embeddings_path)
            embeddings_path="$2"
            shift
            shift
            ;;
        --train_embeddings_path)
            train_embeddings_path="$2"
            shift
            shift
            ;;
        --train)
            train_mode=true
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


if [[ -z "$embeddings_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

if $train_mode && [[ -n "$train_embeddings_path" ]]; then
    echo "Error: --train and --train_embeddings_path cannot be used together."
    exit 1
fi

# Convert embeddings_path to absolute path if it's relative
if [[ "$embeddings_path" != /* ]]; then
    embeddings_path=$(cd "$(dirname "$embeddings_path")" && pwd)/$(basename "$embeddings_path")
fi

# Convert train_embeddings_path to absolute path if it's set and relative
if [[ -n "$train_embeddings_path" ]] && [[ "$train_embeddings_path" != /* ]]; then
    train_embeddings_path=$(cd "$(dirname "$train_embeddings_path")" && pwd)/$(basename "$train_embeddings_path")
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/../src/sequence_metrics"

if [[ -n "$train_embeddings_path" ]] && [[ "$train_embeddings_path" != "" ]]; then
    julia run_min_embedding_distance.jl "$embeddings_path" "$train_embeddings_path"
else
    julia run_min_embedding_distance.jl "$embeddings_path"
fi

if $train_mode; then
    # Copy self results to non-self results
    cp "${embeddings_path%.csv}_min_embedding_distance_self.csv" "${embeddings_path%.csv}_min_embedding_distance.csv"
    echo "Copied self results to non-self results: ${embeddings_path%.csv}_min_embedding_distance_self.csv -> ${embeddings_path%.csv}_min_embedding_distance.csv"
fi
