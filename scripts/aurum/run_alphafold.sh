#!/bin/bash

USAGE="--working_directory <working_directory> --sequence_id <sequence_id> --sequence <sequence>"

############################################################
# Argument parsing                                         #
############################################################
Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --working_directory         Working directory for outputs"
    echo "  --sequence_id               Identifier for the sequence"
    echo "  --sequence                  Amino acid sequence"
    echo "  -h, --help                  Show this help message and exit"
    echo
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --working_directory)
            WRK_DIR="$2"
            shift 2
            ;;
        --sequence_id)
            SEQUENCE_ID="$2"
            shift 2
            ;;
        --sequence)
            SEQUENCE="$2"
            shift 2
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

if [[ -z "$WRK_DIR" ]] || [[ -z "$SEQUENCE_ID" ]] || [[ -z "$SEQUENCE" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

SCRIPT_DIR=$(dirname "$0")
JOBS_DIR="$SCRIPT_DIR"/jobs

LOGS_DIR="$SCRIPT_DIR"/../../logs
mkdir -p "$LOGS_DIR"
OUTPUT_ARG="--output=$LOGS_DIR/%x.%j.out"

sbatch "$OUTPUT_ARG" "$JOBS_DIR"/alphafold.sh "$WRK_DIR" "$SEQUENCE_ID" "$SEQUENCE"
