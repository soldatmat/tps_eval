#!/bin/bash
#
# Matous Soldat, 2026
#
# Compute mean-pooled SaProt structure-aware per-protein embeddings for a directory
# of PDB structures. SA tokens (AA + foldseek 3Di) are produced by SaProt's
# get_struc_seq utility, which shells out to the foldseek binary.
#
# This is a standalone runner (NOT wired into run_eval_pipeline.py). The SaProt env,
# foldseek binary, and SaProt-repo paths default to the Karolina shared-project
# install but can be overridden via flags / env vars.

USAGE="--structs_dir <dir> --output_csv <csv> [--ids_csv <csv>] [--foldseek <bin>] [--saprot_repo <dir>] [--saprot_env <prefix>] [--nogpu]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir   Directory of PDB structures named <ID>.pdb (required)"
    echo "  --output_csv    Output CSV path (required)"
    echo "  --ids_csv       Optional CSV with an Enzyme_marts_ID column to restrict + score coverage"
    echo "  --foldseek      foldseek binary (default: \$FOLDSEEK_PATH or Karolina dplm env)"
    echo "  --saprot_repo   cloned SaProt repo (default: \$SAPROT_REPO or Karolina documents)"
    echo "  --saprot_env    conda env prefix/name (default: \$SAPROT_ENV or Karolina saprot)"
    echo "  --model         HF model id (default: westlake-repl/SaProt_650M_AF2)"
    echo "  --nogpu         Do not use GPU even if available"
    echo "  -h, --help      Show this help message and exit"
}

# Karolina-install defaults (override with flags or env vars).
SAPROT_ENV="${SAPROT_ENV:-/mnt/proj2/fta-26-15/.conda/envs/saprot}"
SAPROT_REPO="${SAPROT_REPO:-/mnt/proj2/fta-26-15/documents/SaProt}"
FOLDSEEK_PATH="${FOLDSEEK_PATH:-/mnt/proj2/fta-26-15/.conda/envs/dplm/bin/foldseek}"
MODEL="westlake-repl/SaProt_650M_AF2"
NOGPU=""
IDS_ARG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --structs_dir) structs_dir="$2"; shift 2;;
        --output_csv) output_csv="$2"; shift 2;;
        --ids_csv) ids_csv="$2"; shift 2;;
        --foldseek) FOLDSEEK_PATH="$2"; shift 2;;
        --saprot_repo) SAPROT_REPO="$2"; shift 2;;
        --saprot_env) SAPROT_ENV="$2"; shift 2;;
        --model) MODEL="$2"; shift 2;;
        --nogpu) NOGPU="--nogpu"; shift;;
        -h|--help) Help; exit 0;;
        *) echo "Unknown option: $1"; Help; exit 1;;
    esac
done

if [ -z "$structs_dir" ] || [ -z "$output_csv" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi
if [ -n "$ids_csv" ]; then
    IDS_ARG="--ids_csv $ids_csv"
fi

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."

eval "$(conda shell.bash hook)"
conda activate "$SAPROT_ENV"
# Karolina compute-node libstdc++/GLIBCXX_3.4.29 fix: prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $CONDA_PREFIX"
echo "Using python: $(which python)"
echo "foldseek: $FOLDSEEK_PATH"
echo "SaProt repo: $SAPROT_REPO"

python src/saprot/extract_saprot_embeddings.py \
    --structs_dir "$structs_dir" \
    --output_csv "$output_csv" \
    --foldseek "$FOLDSEEK_PATH" \
    --saprot_repo "$SAPROT_REPO" \
    --model_location "$MODEL" \
    $IDS_ARG \
    $NOGPU
