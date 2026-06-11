#!/bin/bash

# Login-node FAN-OUT driver wrapper for the orchestrator's `--fold alphafold3`.
# Converts the generated FASTA into an (ID,sequence) CSV, then runs the existing
# src/alphafold/run_alphafold_jobs.py, which submits ONE AlphaFold3 SLURM job per
# sequence (skipping designs whose <ID>.pdb already exists) and PRINTS the list of
# submitted job ids as its last line:  AlphaFold job IDs: ['123', '456']
# The orchestrator's Engine runs THIS script directly (not via sbatch) and parses that
# line for the afterok dependencies of the structure branch. Apo (protein-only) fold;
# ligand/ion co-folding is a future pass (run_alphafold_jobs.py supports it via columns).

USAGE="--cluster <c> --fasta_path <fasta> --working_directory <dir> [--cofold none|mg|mg_ppi] [--model_seeds S1 S2 ...]"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --cluster            Cluster name passed through to submit_job.sh (e.g. aurum)"
    echo "  --fasta_path         Generated sequences FASTA to fold (one AF3 job per record)"
    echo "  --working_directory  AF3 work dir; structures land in <dir>/structs, AF3 trees in <dir>/af_output"
    echo "  --cofold             Co-fold the class-I TPS active site (HOLO): none (default, apo),"
    echo "                       mg (trinuclear Mg2+ cluster), mg_ppi (Mg2+ + diphosphate head group)"
    echo "  --model_seeds        AF3 model seeds (all remaining tokens; default 42)"
    echo "  -h, --help           Show this help message and exit"
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi

model_seeds=()
cofold="none"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cluster) cluster="$2"; shift 2 ;;
        --fasta_path) fasta_path="$2"; shift 2 ;;
        --working_directory) working_directory="$2"; shift 2 ;;
        --cofold) cofold="$2"; shift 2 ;;
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
case "$cofold" in none|mg|mg_ppi) ;; *) echo "Unknown --cofold: $cofold (use none|mg|mg_ppi)"; exit 1 ;; esac
[[ ${#model_seeds[@]} -eq 0 ]] && model_seeds=(42)

# Absolute paths (this runs on the login node from an arbitrary cwd).
[[ "$fasta_path" != /* ]] && fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
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

# FASTA -> (ID,sequence[, co-fold ion columns]) CSV for run_alphafold_jobs.py.
# Co-fold components are class-I TPS active-site ions, all by PDB CCD code (AF3 treats
# ions as ligands w/ ccdCodes; CCD preferred over SMILES for correct geometry):
#   mg     = trinuclear Mg2+ cluster (3x MG, ligated by DDXXD + NSE/DTE).
#   mg_ppi = the cluster + a diphosphate head group (POP, pyrophosphate 2-) that bridges
#            it -- a substrate-agnostic stand-in for the prenyl-PP diphosphate.
# NOTE: AF3 free-ion placement is a HYPOTHESIS, not ground truth -- verify the Mg land at
# the DDXXD/NSE side chains downstream (Christianson 2017; AF3 cofolding caveats).
csv_path="$working_directory/af3_input.csv"
cols_file="$working_directory/.af3_cofold_cols.sh"
python - "$fasta_path" "$csv_path" "$cofold" "$cols_file" <<'PY'
import csv, sys
fa, out, cofold, cols_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
# (ion_id, ccd_code) pairs per preset -> one AF3 ion (ligand w/ ccdCodes) entry each.
COFOLD_IONS = {
    "none":   [],
    "mg":     [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG")],
    "mg_ppi": [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG"), ("PPI", "POP")],
}
ions = COFOLD_IONS[cofold]
id_cols = [f"ion{i+1}_id" for i in range(len(ions))]
ccd_cols = [f"ion{i+1}_ccd" for i in range(len(ions))]
recs, cur, seq = [], None, []
with open(fa) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur is not None:
                recs.append((cur, "".join(seq)))
            cur, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
if cur is not None:
    recs.append((cur, "".join(seq)))
header = ["ID", "sequence"]
for idc, cdc in zip(id_cols, ccd_cols):
    header += [idc, cdc]
ion_vals = [v for (iid, ccd) in ions for v in (iid, ccd)]
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    for rid, s in recs:
        w.writerow([rid, s] + ion_vals)
# Hand the ion column-name lists back to bash (single source of truth = this script).
with open(cols_file, "w") as fh:
    fh.write(f'ION_ID_COLS="{" ".join(id_cols)}"\n')
    fh.write(f'ION_CCD_COLS="{" ".join(ccd_cols)}"\n')
print(f"[fanout] wrote {len(recs)} sequence(s) to {out} (cofold={cofold}, "
      f"{len(ions)} ion(s)/design)")
PY

# Build the ion column args for run_alphafold_jobs.py (empty when --cofold none).
. "$cols_file"
ion_args=()
if [[ -n "${ION_ID_COLS:-}" ]]; then
    ion_args+=(--ion_id_column_names $ION_ID_COLS --ion_ccdcodes_column_names $ION_CCD_COLS)
fi

# Fan out: one AF3 job per sequence. --use_protein_id_as_filename -> <ID>.pdb (filenames
# stay the sequence id regardless of co-fold). run_alphafold_jobs.py prints
# "AlphaFold job IDs: [...]" as its LAST line (the orchestrator parses it).
python src/alphafold/run_alphafold_jobs.py \
    --csv_path "$csv_path" \
    --working_directory "$working_directory" \
    --use_protein_id_as_filename \
    --cluster "$cluster" \
    --model_seeds "${model_seeds[@]}" \
    "${ion_args[@]}"
