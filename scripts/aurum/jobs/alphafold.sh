#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=8
#SBATCH -p b32_128_gpu
#SBATCH --constraint=alphafold3
#SBATCH --mem=100G
#SBATCH --gres=gpu:1

# Usage: sbatch alphafold.sh <working_directory> <proteins> <ligands> [<save_directory>]

############################################################
# Argument parsing                                         #
############################################################
USAGE="--working_directory <working_directory> --sequence_id <sequence_id> --proteins <ID1 SEQ1 ID2 SEQ2 ...> --ligands <ID1 SMILES1 ID2 SMILES2 ...> --ions <ID1 CCDCODE1 ID2 CCDCODE2 ...> --save_directory <save_directory> --model_seeds <SEED1 SEED2 ...>"

Help()
{
    # Display Help
    echo "Usage: alphafold.sh $USAGE"
    echo
    echo "Arguments:"
    echo "  --working_directory         Name of the working directory"
    echo "  --sequence_id               Used to name result files"
    echo "  --proteins                  List of proteins in format: ID1 SEQ1 ID2 SEQ2 ... All following tokens (until next --option) are parsed as proteins."
    echo "  --ligands                   List of ligands in format: ID1 SMILES1 ID2 SMILES2 ... All following tokens (until next --option) are parsed as ligands."
    echo "  --ions                      List of ions in format: ID1 CCDCODE1 ID2 CCDCODE2 ... All following tokens (until next --option) are parsed as ions."
    echo "  --save_directory            Directory to save final pdb structures."
    echo "  --model_seeds               Model seeds to use, separated by space. All following tokens (until next --option) are parsed as model seeds."
    echo "  -h, --help                  Show this help message and exit"
    echo
}

# Parse long options manually
JOB_ARGS=""
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
        --proteins)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                PROTEINS="${PROTEINS:+$PROTEINS }$1"
                shift
            done
            ;;
        --ligands)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                LIGANDS="${LIGANDS:+$LIGANDS }$1"
                shift
            done
            ;;
        --ions)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                IONS="${IONS:+$IONS }$1"
                shift
            done
            ;;
        --save_directory)
            SAVE_DIR="$2"
            shift 2
            ;;
        --model_seeds)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                MODEL_SEEDS="${MODEL_SEEDS:+$MODEL_SEEDS }$1"
                shift
            done
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



############################################################
# Main                                                     #
############################################################
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
# Create config JSON                                       #
############################################################
JSON_FILE="$SEQUENCE_ID.json"
JSON_PATH="${JSON_DIR}/${JSON_FILE}"

eval "$(conda shell.bash hook)"
conda activate tps_eval

JOB_ARGS="--sequence_id $SEQUENCE_ID --proteins $PROTEINS --save_path $JSON_PATH --model_seeds $MODEL_SEEDS"
if [ -n "$LIGANDS" ]; then
    JOB_ARGS+=" --ligands $LIGANDS"
fi
if [ -n "$IONS" ]; then
    JOB_ARGS+=" --ions $IONS"
fi
python ../../../src/alphafold/prepare_input.py $JOB_ARGS


############################################################
# Run AlphaFold3                                           #
############################################################
# alphafold 3 installation folder. Only available in b032.
AF3_DIR="/hpcg/local/soft/alphafold3/"
AF3_SIF="alphafold3-20250108.sif"

echo "Running AlphaFold3 for sequence ${SEQUENCE_ID}"
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
python ../../../vendor/cif_to_pdb/cif_to_pdb.py \
    --input_cif "$STRUCT_PATH" \
    --output_pdb "$STRUCT_SAVE_PATH"
