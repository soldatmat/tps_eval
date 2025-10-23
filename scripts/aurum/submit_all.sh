#!/bin/bash

############################################################
# Argument parsing                                         #
############################################################
USAGE="--fasta_path <fasta_path> [--train_path <train_path> --train_embeddings_path <train_embeddings_path> --structs_dir <structs_dir> --train_structs_dir <train_structs_dir>]"

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
. "../paths.sh" # Load TPS_EVAL_ROOT, TPS_EVAL_ENV variables

SCRIPTS_DIR="$TPS_EVAL_ROOT/scripts"
SRC_DIR="$TPS_EVAL_ROOT/src"
JOBS_DIR="$SCRIPTS_DIR/aurum/jobs"

LOGS_DIR="$TPS_EVAL_ROOT/logs/submit_all-$(date '+%Y_%m_%d_%H%M%S')"
if [[ -d "$LOGS_DIR" ]]; then
    i=1
    while [[ -d "${LOGS_DIR}_$i" ]]; do
        ((i++))
    done
    LOGS_DIR="${LOGS_DIR}_$i"
fi
mkdir -p "$LOGS_DIR"
OUTPUT_ARG="--output=$LOGS_DIR/%x.%j.out"

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"



############################################################
# Prepare CSV                                              #
############################################################
csv_path="${fasta_path%.fasta}.csv"
if [[ -f "$csv_path" ]]; then
    echo "CSV file already exists: $csv_path"
else
    python "$SRC_DIR/prepare_csv.py" --fasta_path "$fasta_path" --csv_path "$csv_path"
fi



############################################################
# AlphaFold structures                                     #
############################################################
# generated data
if [[ -n "$structs_dir" ]] && [[ "$structs_dir" != "" ]]; then
    structs_working_directory="$(dirname "$structs_dir")"

    run_alphafold_jobs_ret=$( \
    python "$SCRIPTS_DIR/run_alphafold_jobs.py" \
        --csv_path "$csv_path" \
        --working_directory "$structs_working_directory" \
        --save_directory "$structs_dir" \
        --submit_args "\"$OUTPUT_ARG\"" \
    )
    echo "$run_alphafold_jobs_ret"

    # Extract job ids from the last [...] in the output and turn into a bash array
    bracket_content=$(echo "$run_alphafold_jobs_ret" | tr '\n' ' ' | sed -n 's/.*\[\([^]]*\)\].*/\1/p')
    if [[ -n "$bracket_content" ]]; then
        # remove whitespace, remove single/double quotes, replace commas with spaces
        alphafold_job_ids=$(echo "$bracket_content" | sed "s/[[:space:]'\" ]//g" | tr ',' ' ')
    fi
    # echo "AlphaFold job IDs: ${alphafold_job_ids[*]}" # debug
fi

# train data
# TODO: check existence of pdb files, try download, run alphafold jobs



############################################################
# Motifs                                                   #
############################################################
# generated data
motifs_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_motifs.csv"
if [[ -f "$motifs_path" ]]; then
    echo "Motifs file already exists: $motifs_path"
else
    echo "Submitting motif search..."
    motif_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/motif_search.sh --fasta_path "$fasta_path")
    echo "$motif_sbatch_ret"
    motif_job_id=${motif_sbatch_ret##* }
fi

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_motifs_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_motifs.csv"
    if [[ -f "$train_motifs_path" ]]; then
        echo "Train motifs file already exists: $train_motifs_path"
    else
        echo "Submitting motif search for train data..."
        train_motif_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/motif_search.sh --fasta_path "$train_path")
        echo "$train_motif_sbatch_ret"
        train_motif_job_id=${train_motif_sbatch_ret##* }
    fi
fi



############################################################
# Min embedding distance to train data                     #
############################################################
embeddings_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_embedding_esm1b.csv"
if [[ -f "$embeddings_path" ]]; then
    echo "Embeddings file already exists: $embeddings_path"
else
    echo "Submitting ESM-1b embedding extraction..."
    esm_embedding_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/esm_embedding.sh --fasta_path "$fasta_path")
    echo "$esm_embedding_sbatch_ret"
    esm_embedding_job_id=${esm_embedding_sbatch_ret##* }
fi

if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    # If `train_embeddings_path` is not provided, check if the embeddings file exists in the same path as `train_path`.
    if [[ -z "$train_embeddings_path" ]] || [[ "$train_embeddings_path" == "" ]]; then
        potential_train_embeddings_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_embedding_esm1b.csv"
        if [[ -f "$potential_train_embeddings_path" ]]; then
            train_embeddings_path="$potential_train_embeddings_path"
            echo "Train embeddings file already exists: $train_embeddings_path"
        fi
    fi

    # If `train_embeddings_path` is still not set, generate train embeddings.
    if [[ -z "$train_embeddings_path" ]] || [[ "$train_embeddings_path" == "" ]]; then
        echo "Submitting ESM-1b embedding extraction for train data..."
        train_esm_embedding_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/esm_embedding.sh --fasta_path "$train_path")
        echo "$train_esm_embedding_sbatch_ret"
        train_esm_embedding_job_id=${train_esm_embedding_sbatch_ret##* }
        train_embeddings_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_embedding_esm1b.csv"
    fi

    # Prepare min_embedding_distance dependencies
    min_embedding_distance_dependency_args=""
    if [[ -n "$train_esm_embedding_job_id" ]]; then
        min_embedding_distance_dependency_args+="--dependency=afterok:$train_esm_embedding_job_id "
    fi
    if [[ -n "$esm_embedding_job_id" ]]; then
        min_embedding_distance_dependency_args+="--dependency=afterok:$esm_embedding_job_id "
    fi

    # Run min_embedding_distance of generated data to train data
    min_embedding_distance_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_embedding_esm1b_min_embedding_distance.csv"
    if [[ -f "$min_embedding_distance_path" ]]; then
        echo "Min embedding distance file already exists: $min_embedding_distance_path"
    else
        echo "Submitting min embedding distance job..."
        min_embedding_distance_sbatch_ret=$(\
        sbatch "$OUTPUT_ARG" \
            $min_embedding_distance_dependency_args \
            "$JOBS_DIR"/min_embedding_distance.sh \
            --embeddings_path "$embeddings_path" \
            --train_embeddings_path "$train_embeddings_path" \
        )
        echo "$min_embedding_distance_sbatch_ret"
        min_embedding_distance_job_id=${min_embedding_distance_sbatch_ret##* }
    fi
fi



############################################################
# Max sequence identity to train data                      #
############################################################
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    max_sequence_identity_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_max_sequence_identity.csv"
    if [[ -f "$max_sequence_identity_path" ]]; then
        echo "Max sequence identity file already exists: $max_sequence_identity_path"
    else
        echo "Submitting max sequence identity job..."
        max_sequence_identity_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh \
                --fasta_path "$fasta_path" \
                --train_path "$train_path" \
            )
        echo "$max_sequence_identity_sbatch_ret"
        max_sequence_identity_job_id=${max_sequence_identity_sbatch_ret##* }
    fi
fi



############################################################
# Min embedding distance self                              #
############################################################
# generated data
min_embedding_distance_self_dependency_args=""
if [[ -n "$esm_embedding_job_id" ]]; then
    min_embedding_distance_self_dependency_args+="--dependency=afterok:$esm_embedding_job_id "
fi

min_embedding_distance_self_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_embedding_esm1b_min_embedding_distance_self.csv"
if [[ -f "$min_embedding_distance_self_path" ]]; then
    echo "Min embedding distance self file already exists: $min_embedding_distance_self_path"
else
    echo "Submitting min embedding distance self job..."
    min_embedding_distance_self_sbatch_ret=$(\
        sbatch "$OUTPUT_ARG" \
            $min_embedding_distance_self_dependency_args \
            "$JOBS_DIR"/min_embedding_distance.sh --embeddings_path "$embeddings_path" \
    )
    echo "$min_embedding_distance_self_sbatch_ret"
    min_embedding_distance_self_job_id=${min_embedding_distance_self_sbatch_ret##* }
fi

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
        echo "Submitting min embedding distance self job for train data..."
        train_min_embedding_distance_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" \
                $train_embedding_min_distance_dependency_args \
                "$JOBS_DIR"/min_embedding_distance.sh --embeddings_path "$train_embeddings_path" --train \
        )
        echo "$train_min_embedding_distance_sbatch_ret"
        train_min_embedding_distance_job_id=${train_min_embedding_distance_sbatch_ret##* }
    fi
fi



############################################################
# Max sequence identity self                               #
############################################################
# generated data
max_sequence_identity_self_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_max_sequence_identity_self.csv"
if [[ -f "$max_sequence_identity_self_path" ]]; then
    echo "Max sequence identity self file already exists: $max_sequence_identity_self_path"
else
    echo "Submitting max sequence identity self job..."
    max_sequence_identity_self_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh --fasta_path "$fasta_path")
    echo "$max_sequence_identity_self_sbatch_ret"
    max_sequence_identity_self_job_id=${max_sequence_identity_self_sbatch_ret##* }
fi

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_max_sequence_identity_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_max_sequence_identity_self.csv"
    if [[ -f "$train_max_sequence_identity_path" ]]; then
        echo "Train max sequence identity file already exists: $train_max_sequence_identity_path"
    else
        echo "Submitting max sequence identity self job for train data..."
        train_max_sequence_identity_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/max_sequence_identity.sh --fasta_path "$train_path" --train)
        echo "$train_max_sequence_identity_sbatch_ret"
        train_max_sequence_identity_job_id=${train_max_sequence_identity_sbatch_ret##* }
    fi
fi



############################################################
# Soluprot                                                 #
############################################################
# generated data
soluprot_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_soluprot.csv"
if [[ -f "$soluprot_path" ]]; then
    echo "Soluprot file already exists: $soluprot_path"
else
    echo "Submitting Soluprot job..."
    soluprot_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/soluprot.sh --fasta_path "$fasta_path")
    echo "$soluprot_sbatch_ret"
    soluprot_job_id=${soluprot_sbatch_ret##* }
fi

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_soluprot_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_soluprot.csv"
    if [[ -f "$train_soluprot_path" ]]; then
        echo "Train Soluprot file already exists: $train_soluprot_path"
    else
        echo "Submitting Soluprot job for train data..."
        train_soluprot_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/soluprot.sh --fasta_path "$train_path")
        echo "$train_soluprot_sbatch_ret"
        train_soluprot_job_id=${train_soluprot_sbatch_ret##* }
    fi
fi



############################################################
# EnzymeExplorer sequence only                             #
############################################################
# generated data
enzyme_explorer_sequence_only_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_enzyme_explorer_sequence_only.csv"
if [[ -f "$enzyme_explorer_sequence_only_path" ]]; then
    echo "EnzymeExplorer sequence only file already exists: $enzyme_explorer_sequence_only_path"
else
    echo "Submitting EnzymeExplorer sequence only job..."
    enzyme_explorer_sequence_only_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer_sequence_only.sh --fasta_path "$fasta_path")
    echo "$enzyme_explorer_sequence_only_sbatch_ret"
    enzyme_explorer_sequence_only_job_id=${enzyme_explorer_sequence_only_sbatch_ret##* }
fi

# train data
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_enzyme_explorer_sequence_only_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_enzyme_explorer_sequence_only.csv"
    if [[ -f "$train_enzyme_explorer_sequence_only_path" ]]; then
        echo "Train EnzymeExplorer sequence only file already exists: $train_enzyme_explorer_sequence_only_path"
    else
        echo "Submitting EnzymeExplorer sequence only job for train data..."
        train_enzyme_explorer_sequence_only_sbatch_ret=$(sbatch "$OUTPUT_ARG" "$JOBS_DIR"/enzyme_explorer_sequence_only.sh --fasta_path "$train_path")
        echo "$train_enzyme_explorer_sequence_only_sbatch_ret"
        train_enzyme_explorer_sequence_only_job_id=${train_enzyme_explorer_sequence_only_sbatch_ret##* }
    fi
fi



############################################################
# EnzymeExplorer (with structures)                         #
############################################################
# generated data
enzyme_explorer_dependency_args=""
if [[ -n "$alphafold_job_ids" ]]; then
    for af_id in $alphafold_job_ids; do
        enzyme_explorer_dependency_args+="--dependency=afterok:$af_id "
    done
fi

enzyme_explorer_struct=false
if [[ -n "$structs_dir" ]] && [[ "$structs_dir" != "" ]]; then
    enzyme_explorer_struct=true
fi
if $enzyme_explorer_struct; then
    enzyme_explorer_path="$(dirname "$fasta_path")/$(basename "$fasta_path" .fasta)_enzyme_explorer.csv"
    if [[ -f "$enzyme_explorer_path" ]]; then
        echo "EnzymeExplorer (with structures) file already exists: $enzyme_explorer_path"
    else
        echo "Submitting EnzymeExplorer (with structures) job..."
        enzyme_explorer_struct_sbatch_ret=$(\
            sbatch "$OUTPUT_ARG" \
                $enzyme_explorer_dependency_args \
                "$JOBS_DIR"/enzyme_explorer.sh \
                --fasta_path "$fasta_path" \
                --structs_dir "$structs_dir" \
            )
        echo "$enzyme_explorer_struct_sbatch_ret"
        enzyme_explorer_struct_job_id=${enzyme_explorer_struct_sbatch_ret##* }
    fi
fi

# train data
train_enzyme_explorer_struct=false
if [[ -n "$train_path" ]] && [[ "$train_path" != "" ]]; then
    train_enzyme_explorer_struct_path="$(dirname "$train_path")/$(basename "$train_path" .fasta)_enzyme_explorer.csv"
    if [[ -f "$train_enzyme_explorer_struct_path" ]]; then
        echo "Train EnzymeExplorer (with structures) file already exists: $train_enzyme_explorer_struct_path"
        train_enzyme_explorer_struct=true
    else if [[ -n "$train_structs_dir" ]] && [[ "$train_structs_dir" != "" ]]; then
        echo "Submitting EnzymeExplorer (with structures) job for train data..."
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



############################################################
# Plots                                                    #
############################################################
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

    echo "Submitting plots for generated and train data..."
    sbatch "$OUTPUT_ARG" \
    $dependency_args \
    "$JOBS_DIR"/plots.sh \
    --fasta_paths "$train_path" "$fasta_path" \
    --data_names "train" "generated" \
    --data_colors "dodgerblue3" "goldenrod1" \
    --save_dir "$plot_save_dir"
else
    echo "Submitting plots for target data only..."
    sbatch "$OUTPUT_ARG" \
    $dependency_args \
    "$JOBS_DIR"/plots.sh \
    --fasta_paths "$fasta_path" \
    --data_names "generated" \
    --data_colors "goldenrod1" \
    --save_dir "$plot_save_dir"
fi
