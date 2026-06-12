#!/bin/bash
# Codon-optimize protein designs and add Golden Gate overhangs, producing an order-ready
# CSV + txt. Standalone ordering utility (NOT part of the eval / submit-all pipeline);
# runs on a login node. See src/order_preparation/ for the logic.

USAGE="<input.fasta|input.csv> [-o <output_prefix>] [--organism <name>] [--overhang_type <type>] [--seq-column <col>] [--id-column <col>]   OR   --sequence <AA> [--id <name>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Codon-optimizes each amino-acid design (default organism: yeast / S. cerevisiae),"
    echo "removes internal BsaI/BsmBI sites, caps homopolymers + holds a GC window, then adds"
    echo "the fixed Golden Gate overhangs (default: Type 3, heterologous protein expression)."
    echo "Writes <prefix>_order.csv and <prefix>_order.txt (id,sequence per line)."
    echo
    echo "Arguments:"
    echo "  input_path        FASTA or CSV/TSV of protein designs (required unless --sequence)"
    echo "  --sequence        A single amino-acid sequence instead of an input file"
    echo "  --id              ID for --sequence mode (default: design)"
    echo "  -o, --output_prefix  Output path prefix (default: input without extension)"
    echo "  --organism        Target organism for codon usage (default: yeast)"
    echo "  --overhang_type   Golden Gate overhang type (default: 'Type 3')"
    echo "  --seq-column      Amino-acid column name (CSV input)"
    echo "  --id-column       ID column name (CSV input)"
    echo "  --method          CodonOptimize method (default: match_codon_usage; or use_best_codon)"
    echo "  --max_homopolymer Longest allowed single-nucleotide run (default: 6; 0 disables)"
    echo "  --gc_min/--gc_max GC-window bounds as fractions (default: 0.30 / 0.65)"
    echo "  --gc_window       GC sliding-window size in bp (default: 50; 0 disables)"
    echo "  --seed            RNG seed for reproducible sampling (default: 0; -1 = random)"
    echo "  --max_attempts    Re-optimization tries to clear a Type IIS site before a design"
    echo "                    is marked FAILED + excluded from the .txt (default: 5)"
    echo "  -h, --help        Show this help and exit"
}

# Pass-through argument parsing: collect the positional input and forward the rest.
input_path=""
has_sequence=false
passthrough=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            Help
            exit 0
            ;;
        -o|--output_prefix)
            output_prefix="$2"; shift 2 ;;
        --sequence)
            has_sequence=true; passthrough+=("$1" "$2"); shift 2 ;;
        -*)
            passthrough+=("$1" "$2"); shift 2 ;;
        *)
            input_path="$1"; shift ;;
    esac
done

if [[ -z "$input_path" ]] && ! $has_sequence; then
    echo "Error: provide an input file or --sequence <AA>."
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert input_path to an absolute path if it's given and relative.
if [[ -n "$input_path" ]] && [[ "$input_path" != /* ]]; then
    input_path=$(cd "$(dirname "$input_path")" && pwd)/$(basename "$input_path")
fi
# Same for output_prefix's directory, if given.
if [[ -n "$output_prefix" ]] && [[ "$output_prefix" != /* ]]; then
    output_prefix=$(cd "$(dirname "$output_prefix")" && pwd)/$(basename "$output_prefix")
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Prepend the env's own libstdc++ (matches the other run_*.sh wrappers).
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/order_preparation

out_args=()
if [[ -n "$output_prefix" ]]; then
    out_args=(-o "$output_prefix")
fi

# In --sequence mode there is no positional input_path to forward.
input_args=()
if [[ -n "$input_path" ]]; then
    input_args=("$input_path")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Preparing order from: $input_path"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Preparing order from a single --sequence"
fi

python run_prepare_order.py "${input_args[@]}" "${out_args[@]}" "${passthrough[@]}"
rc=$?
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished order preparation (rc=$rc)."
exit $rc
