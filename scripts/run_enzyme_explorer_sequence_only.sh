#!/bin/bash

USAGE="--fasta_path <fasta_path>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path                Path to the FASTA file (required)"
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


if [[ -z "$fasta_path" ]] && [[ "$fasta_path" != "" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert fasta_path to absolute path if it's relative
fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/.."
. "scripts/paths.sh" # Load ENZYME_EXPLORER_SEQUENCE_ONLY_PATH variable

eval "$(conda shell.bash hook)"
conda activate "$ENZYME_EXPLORER_SEQUENCE_ONLY_ENV"

output_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_enzyme_explorer_sequence_only.csv"

cd "$ENZYME_EXPLORER_SEQUENCE_ONLY_PATH"
python "scripts/easy_predict_sequence_only.py" \
    --input-fasta-path "$fasta_path" \
    --output-csv-path "$output_path" \
    --detection-threshold 0.0 \
    --detect-precursor-synthases
