#!/bin/bash
# Setup for TmProt (protein melting-temperature predictor), used by run_tmprot.sh.
#
# TmProt is a VENDORED submodule (vendor/TmProt): its standalone CLI package
# (vendor/TmProt/tmprot-1.0) is installed EDITABLE into a dedicated conda env.
# The LoRA adapter ships inside the package; the ESM2-650M base model
# (facebook/esm2_t33_650M_UR50D) is downloaded from Hugging Face on first run.
#
# This script:
#   1. creates the TMPROT_ENV conda env from scripts/tmprot_environment.yml
#   2. `pip install -e vendor/TmProt/tmprot-1.0` (installs the `tmprot` command)
#   3. (best-effort) warms the HF cache with the ESM2-650M base model
#
# GOTCHA (same as vendor/aggrescan3d): the editable install points at the
# submodule working tree, so a `git submodule update`/reset on vendor/TmProt
# de-registers it and the `tmprot` command breaks. After any submodule update,
# re-run: conda run -n "$TMPROT_ENV" pip install -e vendor/TmProt/tmprot-1.0
#
# Usage:  ./setup_tmprot.sh
# Assumes `conda` is already available and vendor/TmProt is initialized
# (`git submodule update --init vendor/TmProt`).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TMPROT_ENV="${TMPROT_ENV:-tmprot}"
TMPROT_PACKAGE="$REPO_ROOT/vendor/TmProt/tmprot-1.0"

if [ ! -f "$TMPROT_PACKAGE/setup.py" ]; then
    echo "ERROR: $TMPROT_PACKAGE not found. Initialize the submodule first:"
    echo "  git submodule update --init vendor/TmProt"
    exit 1
fi

echo "[setup_tmprot] 1/3 creating conda env '$TMPROT_ENV' (py3.10 + torch/transformers/peft)"
conda env create -n "$TMPROT_ENV" -f "$SCRIPT_DIR/tmprot_environment.yml"

echo "[setup_tmprot] 2/3 pip install -e $TMPROT_PACKAGE (installs the 'tmprot' CLI)"
conda run -n "$TMPROT_ENV" pip install -e "$TMPROT_PACKAGE"

echo "[setup_tmprot] 3/3 warming HF cache with facebook/esm2_t33_650M_UR50D (best-effort)"
conda run -n "$TMPROT_ENV" python -c "
from transformers import AutoTokenizer, EsmForSequenceClassification
AutoTokenizer.from_pretrained('facebook/esm2_t33_650M_UR50D')
EsmForSequenceClassification.from_pretrained('facebook/esm2_t33_650M_UR50D', num_labels=1)
print('esm2_t33_650M_UR50D cached')
" || echo "  [warn] HF cache warm-up failed (network?); it will download on first run."

cat <<EOF

==============================================================================
TmProt env + editable CLI installed.
  env:      $TMPROT_ENV
  package:  $TMPROT_PACKAGE (editable)

Set in tps_eval/paths.sh:
  TMPROT_ENV="$TMPROT_ENV"

Self-test:
  conda activate $TMPROT_ENV
  printf '>seqA\n%s\n' MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR \\
    > /tmp/tmprot_test.fasta
  tmprot -i /tmp/tmprot_test.fasta -o /tmp/tmprot_out -d ,
  cat /tmp/tmprot_out/tmprot_test.csv

REMEMBER: after any \`git submodule update\` on vendor/TmProt, re-run:
  conda run -n $TMPROT_ENV pip install -e "$TMPROT_PACKAGE"
==============================================================================
EOF
