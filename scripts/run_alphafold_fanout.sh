#!/bin/bash

# Login-node FAN-OUT driver wrapper for the orchestrator's `--fold alphafold3`.
# Converts the generated FASTA into an (ID,sequence) CSV, then runs the existing
# src/alphafold/run_alphafold_jobs.py, which submits ONE AlphaFold3 SLURM job per
# sequence (skipping designs whose <ID>.pdb already exists) and PRINTS the list of
# submitted job ids as its last line:  AlphaFold job IDs: ['123', '456']
# The orchestrator's Engine runs THIS script directly (not via sbatch) and parses that
# line for the afterok dependencies of the structure branch. Apo (protein-only) fold;
# ligand/ion co-folding is a future pass (run_alphafold_jobs.py supports it via columns).

USAGE="--cluster <c> --fasta_path <fasta> --working_directory <dir> [--model_seeds S1 S2 ...]"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --cluster            Cluster name passed through to submit_job.sh (e.g. aurum)"
    echo "  --fasta_path         Generated sequences FASTA to fold (one AF3 job per record)"
    echo "  --working_directory  AF3 work dir; structures land in <dir>/structs, AF3 trees in <dir>/af_output"
    echo "  --model_seeds        AF3 model seeds (all remaining tokens; default 42)"
    echo "  -h, --help           Show this help message and exit"
}

if [[ $# -lt 1 ]]; then Help; exit 1; fi

model_seeds=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cluster) cluster="$2"; shift 2 ;;
        --fasta_path) fasta_path="$2"; shift 2 ;;
        --working_directory) working_directory="$2"; shift 2 ;;
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

# FASTA -> (ID,sequence) CSV for run_alphafold_jobs.py.
csv_path="$working_directory/af3_input.csv"
python - "$fasta_path" "$csv_path" <<'PY'
import csv, sys
fa, out = sys.argv[1], sys.argv[2]
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
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ID", "sequence"])
    w.writerows(recs)
print(f"[fanout] wrote {len(recs)} sequence(s) to {out}")
PY

# Fan out: one AF3 job per sequence. --use_protein_id_as_filename -> <ID>.pdb (no ligand
# suffix). run_alphafold_jobs.py prints "AlphaFold job IDs: [...]" as its LAST line.
python src/alphafold/run_alphafold_jobs.py \
    --csv_path "$csv_path" \
    --working_directory "$working_directory" \
    --use_protein_id_as_filename \
    --cluster "$cluster" \
    --model_seeds "${model_seeds[@]}"
