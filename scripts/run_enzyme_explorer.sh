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
    echo "  --store_intermediate_results (optional, flag, false by default, set to true with --store_intermediate_results)"
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
        --store_intermediate_results)
            store_intermediate_results=1
            extra_args+=(--store-intermediate-results)
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
# Fix for Karolina /lib64/libstdc++.so.6 being too old (missing GLIBCXX_3.4.29
# required by env's pandas). Prepend the env's libstdc++ (6.0.34, has the symbol).
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"


# EnzymeExplorer (revision branch) installs `predict_with_structures` as a
# console script (pip install -e .). It takes a FASTA (or CSV) directly — no
# prepare_csv step — plus a structures dir, and writes an OUTPUT DIRECTORY with
# TWO CSVs: predictions_plm_domains.csv (structure-based) and
# predictions_plm_only_fallback.csv (PLM-only fallback for proteins without
# usable domain features). NOTE: this schema differs from the old single
# <base>_enzyme_explorer.csv (which had an isTPS column); downstream consumers
# that expect the old format must be updated.
input_path="${sequences_csv_path:-$fasta_path}"
output_dir="$(dirname "$input_path")/$(basename "$input_path" | sed 's/\.[^.]*$//')_enzyme_explorer"

echo "Running EnzymeExplorer with structures (predict_with_structures) with the following parameters:"
echo "  sequences:      $input_path"
echo "  structures dir: $structs_dir"
echo "  output dir:     $output_dir"

# Run from the repo dir so the default reference-domains / model bundles under
# data/ resolve.
cd "$ENZYME_EXPLORER_PATH"

ee_args=(--sequences "$input_path" --structures-dir "$structs_dir" --output-dir "$output_dir")
[[ -n "$csv_id_column" ]] && ee_args+=(--id-column "$csv_id_column")
[[ -n "$n_jobs" ]] && ee_args+=(--n-jobs "$n_jobs")
[[ -n "$plm_batch_size" ]] && ee_args+=(--plm-batch-size "$plm_batch_size")

predict_with_structures "${ee_args[@]}"
