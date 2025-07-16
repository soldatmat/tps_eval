#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=8
#SBATCH -p b32_128_gpu
#SBATCH --constraint=alphafold3
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

# Usage: sbatch aplhafold.sh <working_directory> <sequence_id> <sequence>

# User input data
WRK_DIR=$1 # Adjust if needed
SEQUENCE_ID="$2"
SEQUENCE="$3"

echo "Running on $hostname" # Print the node

SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
cd $(dirname "$SCRIPT_PATH")



############################################################
# Prepare folder structure                                 #
############################################################
JSON_DIR="${WRK_DIR}/af_input"
mkdir -p "$JSON_DIR"
OUTPUT_DIR="${WRK_DIR}/af_output"
mkdir -p "$OUTPUT_DIR"



############################################################
# Prepare config JSON                                      #
############################################################
JULIA_ENV=../../../../../terpene_generation/src # TODO load from a config file

JSON_FILE="$SEQUENCE_ID.json"
JSON_PATH="${JSON_DIR}/${JSON_FILE}"
echo "Preparing input JSON for sequence ${SEQUENCE} at ${JSON_PATH}"
julia --project=$JULIA_ENV ../../../src/alphafold/run_prepare_input.jl "$SEQUENCE_ID" "$SEQUENCE" "$JSON_PATH"



############################################################
# Run AlphaFold3                                           #
############################################################
# alphafold 3 installation folder. Only available in b032.
AF3_DIR="/hpcg/local/soft/alphafold3/"
AF3_SIF="alphafold3-20250108.sif"

echo "Running AlphaFold3 for sequence ${SEQUENCE_ID} with sequence ${SEQUENCE}"
time apptainer exec \
     --nv \
     --bind ${WRK_DIR}/af_input:/root/af_input \
     --bind ${WRK_DIR}/af_output:/root/af_output \
     --bind ${AF3_DIR}/models:/root/models \
     --bind ${AF3_DIR}/db:/root/public_databases \
     --bind ${AF3_DIR}/db:/root/public_databases_fallback \
     --bind ${AF3_DIR}/bin/run_alphafold.py:/usr/local/bin/run_alphafold.py \
     ${AF3_DIR}/bin/${AF3_SIF} \
     python /usr/local/bin/run_alphafold.py \
     --json_path=/root/af_input/${JSON_FILE} \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --db_dir=/root/public_databases_fallback \
     --output_dir=/root/af_output



############################################################
# Extract final pdb structure                              #
############################################################
sequence_id_lowercase=$(echo "$SEQUENCE_ID" | tr '[:upper:]' '[:lower:]')
STRUCT_PATH="${WRK_DIR}/af_output/${sequence_id_lowercase}/${sequence_id_lowercase}_model.cif"
STRUCT_SAVE_PATH="$WRK_DIR"/structs/"$SEQUENCE_ID".pdb
echo "Converting CIF to PDB for sequence ${SEQUENCE_ID} from ${$STRUCT_PATH} to ${STRUCT_SAVE_PATH}"
julia --project=$JULIA_ENV ../../../src/alphafold/cif_to_pdb.jl "$STRUCT_PATH"  "$STRUCT_SAVE_PATH"
