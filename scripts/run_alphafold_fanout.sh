#!/bin/bash

# Login-node FAN-OUT driver wrapper for the orchestrator's `--fold alphafold3`.
# Builds the AlphaFold3 input CSV(s) (apo / Mg / Mg+PPi / Mg+substrate / per-design EE
# substrate) via src/alphafold/build_cofold_input.py, then runs the existing
# src/alphafold/run_alphafold_jobs.py for each input group, submitting ONE AF3 SLURM job
# per sequence (skipping designs whose <ID>.pdb already exists). The orchestrator's Engine
# runs THIS script directly (not via sbatch) and parses the SINGLE final line
#   AlphaFold job IDs: [123, 456]
# for the afterok dependencies of the structure branch (it matches the FIRST such line, so
# we suppress the per-group lines and print one combined line at the end).

USAGE="--cluster <c> --fasta_path <fasta> --working_directory <dir> [--cofold MODE] [--enzymeexplorer_csv <csv>] [--model_seeds S1 S2 ...]"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --cluster            Cluster name passed through to submit_job.sh (e.g. aurum)"
    echo "  --fasta_path         Generated sequences FASTA to fold (one AF3 job per record)"
    echo "  --working_directory  AF3 work dir; structures land in <dir>/structs, AF3 trees in <dir>/af_output"
    echo "  --cofold             Co-fold the class-I TPS active site (HOLO):"
    echo "                         none    (default) apo protein only"
    echo "                         mg      trinuclear Mg2+ cluster (3x CCD MG)"
    echo "                         mg_ppi  Mg cluster + bare diphosphate head group (CCD POP)"
    echo "                         mg_gpp|mg_fpp|mg_ggpp|mg_gfpp  Mg + ONE forced prenyl-PP substrate (SMILES), all designs"
    echo "                         mg_ee   Mg + each design's EnzymeExplorer-predicted substrate (needs --enzymeexplorer_csv)"
    echo "  --enzymeexplorer_csv             EnzymeExplorer seq-only CSV (REQUIRED for --cofold mg_ee)"
    echo "  --model_seeds        AF3 model seeds (all remaining tokens; default 42)"
    echo "  -h, --help           Show this help message and exit"
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi

model_seeds=()
cofold="none"
ee_csv=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cluster) cluster="$2"; shift 2 ;;
        --fasta_path) fasta_path="$2"; shift 2 ;;
        --working_directory) working_directory="$2"; shift 2 ;;
        --cofold) cofold="$2"; shift 2 ;;
        --enzymeexplorer_csv) ee_csv="$2"; shift 2 ;;
        --model_seeds)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do model_seeds+=("$1"); shift; done ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

if [[ -z "$cluster" || -z "$fasta_path" || -z "$working_directory" ]]; then
    echo "Missing required argument."; echo "Usage: $0 $USAGE"; exit 1
fi
[[ ${#model_seeds[@]} -eq 0 ]] && model_seeds=(42)

# Absolute paths (this runs on the login node from an arbitrary cwd).
[[ "$fasta_path" != /* ]] && fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
[[ -n "$ee_csv" && "$ee_csv" != /* ]] && ee_csv="$(cd "$(dirname "$ee_csv")" && pwd)/$(basename "$ee_csv")"
[[ "$working_directory" != /* ]] && working_directory="$(mkdir -p "$working_directory" && cd "$working_directory" && pwd)"
mkdir -p "$working_directory"

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"

# Build the AF3 input CSV(s) + manifest (single source of truth for the co-fold encoding).
# AF3 free-ion/ligand placement is a HYPOTHESIS -- verify the Mg/diphosphate land at the
# DDXXD/NSE cage downstream (ion_site_check / substrate_positioning).
build_out=$(python src/alphafold/build_cofold_input.py \
    --fasta "$fasta_path" --cofold "$cofold" --output_dir "$working_directory" \
    ${ee_csv:+--enzymeexplorer_csv "$ee_csv"} 2>&1)
build_rc=$?
echo "$build_out"
if [[ $build_rc -ne 0 ]]; then echo "[fanout] input build failed (cofold=$cofold)"; exit $build_rc; fi
# Pull the fixed column-name vars the builder printed.
eval "$(echo "$build_out" | grep -E '^(ION_ID_COLS|ION_CCD_COLS|LIG_ID_COL|LIG_SMILES_COL)=')"

manifest="$working_directory/af3_cofold_manifest.tsv"

# Fan out: per input group, one AF3 job per sequence. --use_protein_id_as_filename ->
# <ID>.pdb. We capture each group's stdout, strip its own "AlphaFold job IDs:" line (so the
# orchestrator doesn't match a per-group line), and accumulate the ids.
all_ids=()
while IFS=$'\t' read -r csv has_lig n_designs; do
    [[ "$csv" == "csv_path" ]] && continue
    [[ -z "$csv" ]] && continue
    ion_args=()
    if [[ -n "${ION_ID_COLS:-}" ]]; then
        ion_args+=(--ion_id_column_names $ION_ID_COLS --ion_ccdcodes_column_names $ION_CCD_COLS)
    fi
    lig_args=()
    if [[ "$has_lig" == "1" ]]; then
        lig_args+=(--ligand_id_column_names $LIG_ID_COL --ligand_smiles_column_names $LIG_SMILES_COL)
    fi
    echo "[fanout] folding group $(basename "$csv") ($n_designs design(s), ligand=$has_lig)"
    out=$(python src/alphafold/run_alphafold_jobs.py \
        --csv_path "$csv" \
        --working_directory "$working_directory" \
        --use_protein_id_as_filename \
        --cluster "$cluster" \
        --model_seeds "${model_seeds[@]}" \
        "${ion_args[@]}" "${lig_args[@]}")
    echo "$out" | grep -v "AlphaFold job IDs:"
    for id in $(echo "$out" | grep "AlphaFold job IDs:" | grep -oE "[0-9]+"); do
        all_ids+=("$id")
    done
done < "$manifest"

# Single combined contract line for the orchestrator (the FIRST/only "AlphaFold job IDs:").
echo "AlphaFold job IDs: [$(IFS=,; echo "${all_ids[*]}")]"
