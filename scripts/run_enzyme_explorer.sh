#!/bin/bash

USAGE="[--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path            Path to the FASTA file with sequences (required one of: --sequences_csv_path or --fasta_path)"
    echo "  --sequences_csv_path    Path to the CSV file with columns ID and sequence (required one of: --sequences_csv_path or --fasta_path)"
    echo "  --structs_dir           Path to the directory containing structures (required)"
    echo "  -h, --help              Show this help message and exit"
    echo
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --sequences_csv_path)
            sequences_csv_path="$2"
            shift
            shift
            ;;
        --fasta_path)
            fasta_path="$2"
            shift
            shift
            ;;
        --structs_dir)
            structs_dir="$2"
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


if { [[ -z "$sequences_csv_path" ]] && [[ -z "$fasta_path" ]]; } || [[ -z "$structs_dir" ]]; then
    echo "Usage: $0 $USAGE"
    echo
    echo "At least one of --sequences_csv_path or --fasta_path must be provided, and --structs_dir is required."
    exit 1
fi

# Convert sequences_csv_path to absolute path if it's relative
if [[ -n "$sequences_csv_path" ]] && [[ "$sequences_csv_path" != "" ]]; then
    sequences_csv_path="$(cd "$(dirname "$sequences_csv_path")" && pwd)/$(basename "$sequences_csv_path")"
fi

# Convert fasta_patha to absolute path if it's relative
if [[ -n "$fasta_patha" ]] && [[ "$fasta_patha" != "" ]]; then
    fasta_patha="$(cd "$(dirname "$fasta_patha")" && pwd)/$(basename "$fasta_patha")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/.."
. "scripts/paths.sh" # Load ENZYME_EXPLORER_PATH, TPS_EVAL_ROOT variables

eval "$(conda shell.bash hook)"
conda activate terpene_miner


if [[ -z "$sequences_csv_path" ]]; then
    sequences_csv_path="${fasta_path%.fasta}.csv"
    python src/enzyme_explorer/prepare_csv.py --fasta_path "$fasta_path" --csv_path "$sequences_csv_path"
fi

echo "Running enzyme explorer (easy_predict_batching.py) with the following parameters:"
echo "  sequences CSV path: $sequences_csv_path"
echo "  structures directory: $structs_dir"

output_path="$(dirname "$sequences_csv_path")/$(basename "$sequences_csv_path" .csv)_enzyme_explorer.csv"

# The easy_predict-batching.py scripts has to be run in EnzymeExplorer/scripts/ directory
cd "$ENZYME_EXPLORER_PATH/scripts"
python "$TPS_EVAL_ROOT/src/enzyme_explorer/easy_predict-batching.py" \
    --input-directory-with-structures "$structs_dir" \
    --needed-proteins-csv-path "$sequences_csv_path" \
    --csv-id-column ID \
    --n-jobs 20 \
    --is-bfactor-confidence \
    --output-csv-path $output_path \
    --detection-threshold 0 \
    --detect-precursor-synthases \
    --plm-batch-size 20
