#!/bin/bash

USAGE="--fasta_path <fasta_path> [--train_path <train_path> --train_embeddings_path <train_embeddings_path> --structs_dir <structs_dir>]"

############################################################
# Parameters                                               #
############################################################
MOTIFS=("DD..D" "(N|D)D(L|I|V).(S|T)...E")

############################################################
# Argument parsing                                         #
############################################################
Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path                Path to the FASTA file (required)"
    echo "  --train_path                Path to the reference FASTA file (optional)"
    echo "  --train_embeddings_path     Path to the reference embeddings file (optional)"
    echo "  --structs_dir               Directory containing structures (optional)"
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
        --train_path)
            train_path="$2"
            shift
            shift
            ;;
        --train_embeddings_path)
            train_embeddings_path="$2"
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


if [[ -z "$fasta_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
JOBS_DIR="$SCRIPT_DIR"/jobs

LOGS_DIR="$SCRIPT_DIR"/../../logs/submit_all-"$(date '+%Y_%m_%d_%H%M%S')"
mkdir -p "$LOGS_DIR"
OUTPUT_ARG="--output=$LOGS_DIR/%x.%j.out"



sbatch "$OUTPUT_ARG" "$JOBS_DIR"/motif_search.sh --fasta_path "$fasta_path" "${MOTIFS[@]}"

esm_embedding_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/esm_embedding.sh --fasta_path "$fasta_path")
echo "$esm_embedding_sbatch_ret"
esm_embedding_job_id=${esm_embedding_sbatch_ret##* }
embeddings_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_embedding_esm1b.csv"

if [[ -n "$train_embeddings_path" ]] && [[ "$train_embeddings_path" != "" ]]; then
    sbatch "$OUTPUT_ARG" --dependency=afterok:$esm_embedding_job_id \
        "$JOBS_DIR"/min_embedding_distance.sh \
        --embeddings_path "$embeddings_path" \
        --train_embeddings_path "$train_embeddings_path"
fi

if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh \
        --fasta_path "$fasta_path" \
        --train_path "$train_path"
fi

sbatch "$OUTPUT_ARG" --dependency=afterok:$esm_embedding_job_id \
    "$JOBS_DIR"/min_embedding_distance.sh --embeddings_path "$embeddings_path"

sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh --fasta_path "$fasta_path"

sbatch "$OUTPUT_ARG" "$JOBS_DIR"/soluprot.sh --fasta_path "$fasta_path"

sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer_sequence_only.sh --fasta_path "$fasta_path"

if [[ -n "$structs_dir" ]] && [[ "$structs_dir" != "" ]]; then
    sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer.sh \
        --fasta_path "$fasta_path" \
        --structs_dir "$structs_dir"
fi
