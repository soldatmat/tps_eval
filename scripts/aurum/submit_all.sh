#!/bin/bash

USAGE="--fasta_path <fasta_path> [--train_path <train_path> --train_embeddings_path <train_embeddings_path> --structs_dir <structs_dir> --train_structs_dir <train_structs_dir>]"

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
    echo "  --train_embeddings_path     Custom path to the reference embeddings file (optional), otherwise <train_path>_embedding_esm1b.csv will be used or the embeddings will be generated."
    echo "  --structs_dir               Directory containing structures (optional)"
    echo "  --train_structs_dir         Directory containing train structures (optional)"
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
        --train_structs_dir)
            train_structs_dir="$2"
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



########## motifs ##########
# generated data
sbatch "$OUTPUT_ARG" "$JOBS_DIR"/motif_search.sh --fasta_path "$fasta_path"

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_motifs_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_motifs.csv"
    if [[ -f "$train_motifs_path" ]]; then
        echo "Train motifs file already exists: $train_motifs_path"
    else
        train_motif_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/motif_search.sh --fasta_path "$train_path")
        echo "$train_motif_sbatch_ret"
        train_motif_job_id=${train_motif_sbatch_ret##* }
    fi
fi



########## min embedding distance to train data ##########
esm_embedding_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/esm_embedding.sh --fasta_path "$fasta_path")
echo "$esm_embedding_sbatch_ret"
esm_embedding_job_id=${esm_embedding_sbatch_ret##* }
embeddings_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_embedding_esm1b.csv"

if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    # If `train_embeddings_path` is not provided, check if the embeddings file exists in the same path as `train_path`.
    if [[ -z "$train_embeddings_path" ]] || [[ "$train_embeddings_path" == "" ]]; then
        potential_train_embeddings_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_embedding_esm1b.csv"
        if [[ -f "$potential_train_embeddings_path" ]]; then
            train_embeddings_path="$potential_train_embeddings_path"
        fi
    fi

    if [[ -n "$train_embeddings_path" ]] && [[ "$train_embeddings_path" != "" ]]; then
        min_embedding_distance_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" \
                --dependency=afterok:$esm_embedding_job_id \
                "$JOBS_DIR"/min_embedding_distance.sh \
                --embeddings_path "$embeddings_path" \
                --train_embeddings_path "$train_embeddings_path" \
        )
    else
        train_esm_embedding_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/esm_embedding.sh --fasta_path "$train_path")
        echo "$train_esm_embedding_sbatch_ret"
        train_esm_embedding_job_id=${train_esm_embedding_sbatch_ret##* }
        train_embeddings_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_embedding_esm1b.csv"

        min_embedding_distance_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" \
                --dependency=afterok:$esm_embedding_job_id \
                --dependency=afterok:$train_esm_embedding_job_id \
                "$JOBS_DIR"/min_embedding_distance.sh \
                --embeddings_path "$embeddings_path" \
                --train_embeddings_path "$train_embeddings_path" \
        )
    fi
    echo "$min_embedding_distance_sbatch_ret"
    min_embedding_distance_job_id=${min_embedding_distance_sbatch_ret##* }
fi



########## max sequence identity to train data ##########
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    max_sequence_identity_sbatch_ret=$(\
        sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh \
            --fasta_path "$fasta_path" \
            --train_path "$train_path" \
        )
    echo "$max_sequence_identity_sbatch_ret"
    max_sequence_identity_job_id=${max_sequence_identity_sbatch_ret##* }
fi



########## min embedding distance self ##########
# generated data
min_embedding_distance_self_sbatch_ret=$(\
    sbatch "$OUTPUT_ARG" --dependency=afterok:$esm_embedding_job_id \
        "$JOBS_DIR"/min_embedding_distance.sh --embeddings_path "$embeddings_path" \
)
echo "$min_embedding_distance_self_sbatch_ret"
min_embedding_distance_self_job_id=${min_embedding_distance_self_sbatch_ret##* }

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_min_embedding_distance_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_embedding_esm1b_min_embedding_distance_self.csv"
    if [[ -f "$train_min_embedding_distance_path" ]]; then
        echo "Train min embedding distance file already exists: $train_min_embedding_distance_path"
    else
        train_embedding_min_distance_dependency_args=""
        if [[ -n "$train_esm_embedding_job_id" ]]; then
            train_embedding_min_distance_dependency_args+="--dependency=afterok:$train_esm_embedding_job_id "
        fi
        # TODO delete testing echoes
        train_min_embedding_distance_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" \
                $train_embedding_min_distance_dependency_args \
                "$JOBS_DIR"/min_embedding_distance.sh --embeddings_path "$train_embeddings_path" --train \
        )
        echo "$train_min_embedding_distance_sbatch_ret"
        train_min_embedding_distance_job_id=${train_min_embedding_distance_sbatch_ret##* }
    fi
fi



########## max sequence identity self ##########
# generated data
max_sequence_identity_self_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh --fasta_path "$fasta_path")
echo "$max_sequence_identity_self_sbatch_ret"
max_sequence_identity_self_job_id=${max_sequence_identity_self_sbatch_ret##* }

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_max_sequence_identity_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_max_sequence_identity_self.csv"
    if [[ -f "$train_max_sequence_identity_path" ]]; then
        echo "Train max sequence identity file already exists: $train_max_sequence_identity_path"
    else
        train_max_sequence_identity_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh --fasta_path "$train_path" --train)
        echo "$train_max_sequence_identity_sbatch_ret"
        train_max_sequence_identity_job_id=${train_max_sequence_identity_sbatch_ret##* }
    fi
fi



########## Soluprot ##########
# generated data
soluprot_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/soluprot.sh --fasta_path "$fasta_path")
echo "$soluprot_sbatch_ret"
soluprot_job_id=${soluprot_sbatch_ret##* }

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_soluprot_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_soluprot.csv"
    if [[ -f "$train_soluprot_path" ]]; then
        echo "Train Soluprot file already exists: $train_soluprot_path"
    else
        train_soluprot_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/soluprot.sh --fasta_path "$train_path")
        echo "$train_soluprot_sbatch_ret"
        train_soluprot_job_id=${train_soluprot_sbatch_ret##* }
    fi
fi



########## Enzyme Explorer sequence only ##########
# generated data
enzyme_explorer_sequence_only_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer_sequence_only.sh --fasta_path "$fasta_path")
echo "$enzyme_explorer_sequence_only_sbatch_ret"
enzyme_explorer_sequence_only_job_id=${enzyme_explorer_sequence_only_sbatch_ret##* }

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_enzyme_explorer_sequence_only_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_enzyme_explorer_sequence_only.csv"
    if [[ -f "$train_enzyme_explorer_sequence_only_path" ]]; then
        echo "Train EnzymeExplorer sequence only file already exists: $train_enzyme_explorer_sequence_only_path"
    else
        train_enzyme_explorer_sequence_only_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer_sequence_only.sh --fasta_path "$train_path")
        echo "$train_enzyme_explorer_sequence_only_sbatch_ret"
        train_enzyme_explorer_sequence_only_job_id=${train_enzyme_explorer_sequence_only_sbatch_ret##* }
    fi
fi


########## Enzyme Explorer (with structures) ##########
# generated data
enzyme_explorer_struct=false
if [[ -n "$structs_dir" ]] && [[ "$structs_dir" != "" ]]; then
    enzyme_explorer_struct=true
fi
if $enzyme_explorer_struct; then
    enzyme_explorer_struct_sbatch_ret=$(\
        sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer.sh \
            --fasta_path "$fasta_path" \
            --structs_dir "$structs_dir" \
        )
    echo "$enzyme_explorer_struct_sbatch_ret"
    enzyme_explorer_struct_job_id=${enzyme_explorer_struct_sbatch_ret##* }
fi

# train data
train_enzyme_explorer_struct=false
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_enzyme_explorer_struct_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_enzyme_explorer.csv"
    if [[ -f "$train_enzyme_explorer_struct_path" ]]; then
        echo "Train EnzymeExplorer (with structures) file already exists: $train_enzyme_explorer_struct_path"
        train_enzyme_explorer_struct=true
    else if [[ -n "$train_structs_dir" ]] && [[ "$train_structs_dir" != "" ]]; then
        train_enzyme_explorer_struct_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer.sh \
                --fasta_path "$train_path" \
                --structs_dir "$train_structs_dir" \
            )
        echo "$train_enzyme_explorer_struct_sbatch_ret"
        train_enzyme_explorer_struct_job_id=${train_enzyme_explorer_struct_sbatch_ret##* }
        train_enzyme_explorer_struct=true
    else
        echo "Train EnzymeExplorer (with structures) file does not exist (in $train_enzyme_explorer_struct_path) and no train_structs_dir provided."
    fi
fi
fi



########## Plots ##########
plot_save_dir=$(dirname "$fasta_path")/plots

# Build dependency string dynamically
dependency_args=""
if [[ -n "$min_embedding_distance_job_id" ]]; then
    dependency_args+="--dependency=afterok:$min_embedding_distance_job_id "
fi
if [[ -n "$min_embedding_distance_self_job_id" ]]; then
    dependency_args+="--dependency=afterok:$min_embedding_distance_self_job_id "
fi
if [[ -n "$max_sequence_identity_job_id" ]]; then
    dependency_args+="--dependency=afterok:$max_sequence_identity_job_id "
fi
if [[ -n "$max_sequence_identity_self_job_id" ]]; then
    dependency_args+="--dependency=afterok:$max_sequence_identity_self_job_id "
fi
if [[ -n "$soluprot_job_id" ]]; then
    dependency_args+="--dependency=afterok:$soluprot_job_id "
fi
if [[ -n "$enzyme_explorer_sequence_only_job_id" ]]; then
    dependency_args+="--dependency=afterok:$enzyme_explorer_sequence_only_job_id "
fi
if [[ -n "$enzyme_explorer_struct_job_id" ]]; then
    dependency_args+="--dependency=afterok:$enzyme_explorer_struct_job_id "
fi

if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    if [[ -n "$train_min_embedding_distance_job_id" ]]; then
        dependency_args+="--dependency=afterok:$train_min_embedding_distance_job_id "
    fi
    if [[ -n "$train_max_sequence_identity_job_id" ]]; then
        dependency_args+="--dependency=afterok:$train_max_sequence_identity_job_id "
    fi
    if [[ -n "$train_soluprot_job_id" ]]; then
        dependency_args+="--dependency=afterok:$train_soluprot_job_id "
    fi
    if [[ -n "$train_enzyme_explorer_sequence_only_job_id" ]]; then
        dependency_args+="--dependency=afterok:$train_enzyme_explorer_sequence_only_job_id "
    fi
    if [[ -n "$train_enzyme_explorer_struct_job_id" ]]; then
        dependency_args+="--dependency=afterok:$train_enzyme_explorer_struct_job_id "
    fi

    sbatch "$OUTPUT_ARG" \
    $dependency_args \
    "$JOBS_DIR"/plots.sh \
    --fasta_paths "$train_path" "$fasta_path" \
    --data_names "train" "generated" \
    --data_colors "dodgerblue3" "goldenrod1" \
    --save_dir "$plot_save_dir"
else
    sbatch "$OUTPUT_ARG" \
    $dependency_args \
    "$JOBS_DIR"/plots.sh \
    --fasta_paths "$fasta_path" \
    --data_names "generated" \
    --data_colors "goldenrod1" \
    --save_dir "$plot_save_dir"
fi
