#!/bin/bash

USAGE="[--sequences_csv_path <sequences_csv_path> --fasta_path <fasta_path>] --structs_dir <structs_dir>"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path                 Path to the FASTA file with sequences (required one of: --sequences_csv_path or --fasta_path)"
    echo "  --sequences_csv_path         Path to the CSV file with columns ID and sequence (required one of: --sequences_csv_path or --fasta_path)"
    echo "  --structs_dir                Path to the directory containing structures (required)"
    echo "  --csv_id_column              Column name for sequence IDs in the CSV file (optional)"
    echo "  --n_jobs                     Number of parallel jobs for prediction (optional)"
    echo "  --is_bfactor_confidence      (optional, flag, true by default, set to false with --no-is_bfactor_confidence)"
    echo "  --detection_threshold        (optional)"
    echo "  --detect_precursor_synthases (optional, flag, true by default, set to false with --no-detect_precursor_synthases)"
    echo "  --plm_batch_size             (optional)"
    echo "  -h, --help                   Show this help message and exit"
    echo
}

# Collect extra arguments for easy_predict.py
extra_args=()

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
        --csv_id_column)
            csv_id_column="$2"
            extra_args+=(--csv-id-column "$csv_id_column")
            shift
            shift
            ;;
        --n_jobs)
            n_jobs="$2"
            extra_args+=(--n-jobs "$n_jobs")
            shift
            shift
            ;;
        --is_bfactor_confidence)
            is_bfactor_confidence=1
            shift
            ;;
        --no-is_bfactor_confidence)
            is_bfactor_confidence=0
            shift
            ;;
        --detection_threshold)
            detection_threshold="$2"
            extra_args+=(--detection-threshold "$detection_threshold")
            shift
            shift
            ;;
        --detect_precursor_synthases)
            detect_precursor_synthases=1
            shift
            ;;
        --no-detect_precursor_synthases)
            detect_precursor_synthases=0
            shift
            ;;
        --plm_batch_size)
            plm_batch_size="$2"
            extra_args+=(--plm-batch-size "$plm_batch_size")
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
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load ENZYME_EXPLORER_ENV, ENZYME_EXPLORER_PATH variables

eval "$(conda shell.bash hook)"
conda activate "$ENZYME_EXPLORER_ENV"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"


if [[ -z "$sequences_csv_path" ]]; then
    sequences_csv_path="${fasta_path%.fasta}.csv"
    python src/enzyme_explorer/prepare_csv.py --fasta_path "$fasta_path" --csv_path "$sequences_csv_path"
fi

echo "Running enzyme explorer (easy_predict.py) with the following parameters:"
echo "  sequences CSV path: $sequences_csv_path"
echo "  structures directory: $structs_dir"

output_path="$(dirname "$sequences_csv_path")/$(basename "$sequences_csv_path" .csv)_enzyme_explorer.csv"

# The easy_predict.py scripts has to be run in EnzymeExplorer/scripts/ directory
cd "$ENZYME_EXPLORER_PATH/scripts"

# Set defaults for these options if not provided by the user
if [[ -z "$csv_id_column" ]]; then
    extra_args+=(--csv-id-column "ID")
fi
if [[ -z "$n_jobs" ]]; then
    extra_args+=(--n-jobs "20")
fi
if [[ -z "$plm_batch_size" ]]; then
    extra_args+=(--plm-batch-size "20")
fi
# Ensure is_bfactor_confidence defaults to true (1) when not provided
if [[ -z "$is_bfactor_confidence" ]]; then
    is_bfactor_confidence=1
fi
if [[ "$is_bfactor_confidence" == "1" ]]; then
    extra_args+=(--is-bfactor-confidence)
fi
# Ensure detect_precursor_synthases defaults to true (1) when not provided
if [[ -z "$detect_precursor_synthases" ]]; then
    detect_precursor_synthases=1
fi
if [[ "$detect_precursor_synthases" == "1" ]]; then
    extra_args+=(--detect-precursor-synthases)
fi

python easy_predict.py \
    --input-directory-with-structures "$structs_dir" \
    --needed-proteins-csv-path "$sequences_csv_path" \
    ${extra_args[@]} \
    --output-csv-path "$output_path"
