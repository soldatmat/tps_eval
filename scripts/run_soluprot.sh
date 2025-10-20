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
cd "$SCRIPT_DIR/.."
. "scripts/paths.sh" # Load SOLUPROT_PATH, SOLUPROT_ENV variables

eval "$(conda shell.bash hook)"
conda activate "$SOLUPROT_ENV"

if [[ -z "$SCRATCH" ]] || [[ ! -d "$SCRATCH" ]]; then
    echo "Error: SCRATCH variable is not set to a valid directory."
    exit 1
fi

python "$SOLUPROT_PATH/soluprot.py" \
    --i_fa "$fasta_path" \
    --o_csv "$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_soluprot.csv" \
    --tmp_dir "$SCRATCH"
