#!/bin/bash
# Login-node dispatcher for the selection layer (merge / select). Unlike the metric
# runners these are CPU-only, seconds-fast pandas ops -- run them directly on the login
# node (no SLURM), so run_funnel.py can call them inline between tiers.
#
# Usage:
#   scripts/run_selection.sh merge  --entries <csv|dir|glob> [...] --output merged.csv
#   scripts/run_selection.sh select --merged merged.csv --spec spec.json \
#                                   --output_prefix phaseN [--fasta seed.fasta]

USAGE="<merge|select> [op args...]"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

op="$1"; shift
case "$op" in
    merge)  entry="run_merge.py"  ;;
    select) entry="run_select.py" ;;
    -h|--help)
        echo "Usage: $0 $USAGE"
        exit 0
        ;;
    *)
        echo "Unknown selection op: $op (expected: merge | select)"
        echo "Usage: $0 $USAGE"
        exit 1
        ;;
esac

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Karolina compute-node libstdc++/GLIBCXX_3.4.29 fix (see other runners). Also needed
# for mmseqs (diversity_dedup) which links the env's libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd src/selection
python "$entry" "$@"
exit $?
