#!/bin/bash

############################################################
# Parameters                                               #
############################################################
# Class-I terpene synthase metal-binding motifs:
#   DD..D                   : canonical "aspartate-rich" motif (strict).
#   D[DE]..[DE]             : relaxed acidic variant of the above — tolerates the
#                             conservative D->E substitution at the 2nd/5th positions
#                             (DExxD / DDxxE / DExxE), which still supplies a carboxylate
#                             for Mg2+ coordination. Generative models (e.g. DPLM run_41_V)
#                             often emit this variant, so the strict DD..D undercounts a
#                             functionally-present first metal site.
#   [DE][DE]..[DE]          : fully-relaxed acidic variant — also allows D->E at the
#                             1st position (EDxxD, EExxE, etc.). Superset of the two
#                             above; the three are nested (strict < semi < fully relaxed)
#                             and kept as separate columns to give graded hits.
#   [DN]D..[STG]...E        : NSE/DTE second metal-binding motif — RELAXED as of
#                             2026-08-18 from the old strict (N|D)D(L|I|V).(S|T)...E.
#                             Source: Durairaj et al. 2019, Phytochemistry 158:157-165
#                             (an expansion of Christianson 2006's Chem. Rev. consensus),
#                             independently confirmed optimal by an empirical grid search
#                             over our own MARTS-DB Class I/II sequences the same day
#                             (Class I recall 58.9% -> 80.8%; Class II false-positive rate
#                             stays low at 2.7%, vs strict DDxxD's own 3.8% reference rate).
#                             Drops the (L|I|V) restriction at position 3 (now any residue)
#                             and adds G alongside S/T at the (S/T) slot; length/offsets
#                             unchanged (still a 9-mer, coordinating residues at 0,1,4,8).
MOTIFS=("DD..D" "D[DE]..[DE]" "[DE][DE]..[DE]" "[DN]D..[STG]...E")

############################################################
# Script                                                   #
############################################################
USAGE="--fasta_path <fasta_path> [<motif1> <motif2> ...]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_path            Path to the FASTA file (required)"
    echo "  <motif1> <motif2> ...   Any number of motifs to search for in the sequences"
    echo "  -h, --help              Show this help message and exit"
    echo
}

# Parse long options manually
motifs=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_path)
            fasta_path="$2"
            shift
            shift
            ;;
        -h|--help)
            Help
            exit 0
            ;;
        *)
            motifs+=("$key")
            shift
            ;;
    esac
done


if [[ -z "$fasta_path" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert fasta_path to absolute path if it's relative
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"
fi

if [[ ${#motifs[@]} -eq 0 ]]; then
    motifs=("${MOTIFS[@]}")
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/../.."
# Every tool activates its own conda env, but only the main tps_eval env gets
# `pip install -e .` -- put the in-repo package on PYTHONPATH so `import tps_eval`
# works in ALL of them (absolute, so it survives later `cd`s and child processes).
export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"
. ./paths.sh # Load TPS_EVAL_ENV

eval "$(conda shell.bash hook)"
conda activate "$TPS_EVAL_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"


python -m tps_eval.sequence_metrics.run_motif_search "$fasta_path" "${motifs[@]}"
