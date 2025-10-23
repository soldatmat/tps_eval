#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=8
#SBATCH -p b32_128_gpu
#SBATCH --constraint=alphafold3
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

# Usage: sbatch aplhafold.sh <working_directory> <sequence_id> <sequence> [<save_directory>]

# User input data
WRK_DIR="$1"
SEQUENCE_ID="$2"
SEQUENCE="$3"
SAVE_DIR="$4"

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
JSON_FILE="$SEQUENCE_ID.json"
JSON_PATH="${JSON_DIR}/${JSON_FILE}"

eval "$(conda shell.bash hook)"
conda activate tps_eval

python ../../../src/alphafold/prepare_input.py \
    --sequence_id "$SEQUENCE_ID" \
    --sequence "$SEQUENCE" \
    --save_path "$JSON_PATH"



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
eval "$(conda shell.bash hook)"
conda activate tps_eval

sequence_id_lowercase=$(echo "$SEQUENCE_ID" | tr '[:upper:]' '[:lower:]')
STRUCT_PATH="${WRK_DIR}/af_output/${sequence_id_lowercase}/${sequence_id_lowercase}_model.cif"
if [ -n "$SAVE_DIR" ]; then
    STRUCT_SAVE_PATH="${SAVE_DIR}/${SEQUENCE_ID}.pdb"
else
    STRUCT_SAVE_PATH="${WRK_DIR}/structs/${SEQUENCE_ID}.pdb"
fi
mkdir -p "$(dirname "$STRUCT_SAVE_PATH")"
echo "Converting CIF to PDB for sequence ${SEQUENCE_ID} from ${STRUCT_PATH} to ${STRUCT_SAVE_PATH}"
python ../../../src/alphafold/cif_to_pdb.py \
    --input_cif "$STRUCT_PATH" \
    --output_pdb "$STRUCT_SAVE_PATH"
