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
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. "./paths.sh" # Load ENZYME_EXPLORER_SEQUENCE_ONLY_PATH, ENZYME_EXPLORER_SEQUENCE_ONLY_ENV variables

eval "$(conda shell.bash hook)"
conda activate "$ENZYME_EXPLORER_SEQUENCE_ONLY_ENV"
# Fix for Karolina /lib64/libstdc++.so.6 being too old (missing GLIBCXX_3.4.29
# required by env's pandas). Prepend the env's libstdc++ (6.0.34, has the symbol).
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"


output_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_enzyme_explorer_sequence_only.csv"

echo "Running EnzymeExplorer sequence-only (predict_sequences_only) with the following parameters:"
echo "  sequences FASTA path: $fasta_path"
echo "  output CSV path: $output_path"

# EnzymeExplorer (revision branch) installs `predict_sequences_only` as a console
# script (pip install -e .). Run from the repo dir so its default model bundles
# under data/ (enzyme_explorer_plm_checkpoints.pkl, calibration_fit_summary.csv)
# resolve. Output schema: id, sequence, <class>_score, <class>_p_calibrated.
cd "$ENZYME_EXPLORER_SEQUENCE_ONLY_PATH"
predict_sequences_only --sequences "$fasta_path" --output-csv "$output_path"
