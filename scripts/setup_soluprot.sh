#!/bin/bash
# Optional setup for SoluProt (protein solubility predictor), used by
# run_soluprot.sh. SoluProt is NOT a pip/conda package: it's a standalone
# download + an old py3.7 conda env + two external binaries (USEARCH, TMHMM).
#
# This script automates the parts that can be automated:
#   1. the py3.7 conda env from scripts/soluprot_environment.yml
#   2. the SoluProt standalone (code + gradient-boosting models)
#   3. a 64-bit USEARCH binary (public domain — works on x86-64 compute nodes;
#      the old usearch11 32-bit binary does NOT)
# TMHMM 2.0 is academic-licensed and must be installed manually (step 4).
#
# Usage:  ./setup_soluprot.sh [install_dir]      (default: $HOME/soluprot_install)
# Assumes `conda` is already available (run from a shell where conda is set up).
# A Karolina-adapted variant lives at scripts/karolina/setup_soluprot.sh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="${1:-$HOME/soluprot_install}"
SOLUPROT_ENV="${SOLUPROT_ENV:-soluprot}"
mkdir -p "$INSTALL_DIR"

echo "[setup_soluprot] 1/3 creating conda env '$SOLUPROT_ENV' (py3.7, scikit-learn 0.20.1)"
conda env create -n "$SOLUPROT_ENV" -f "$SCRIPT_DIR/soluprot_environment.yml"

echo "[setup_soluprot] 2/3 downloading SoluProt standalone"
cd "$INSTALL_DIR"
curl -fsSL 'https://loschmidt.chemi.muni.cz/soluprot/?f=soluprot.zip' -o soluprot.zip
unzip -qo soluprot.zip
rm -f soluprot.zip
SOLUPROT_PATH="$(dirname "$(find "$INSTALL_DIR" -maxdepth 3 -name soluprot.py | head -1)")"
[ -n "$SOLUPROT_PATH" ] || { echo "ERROR: soluprot.py not found after unzip"; exit 1; }

echo "[setup_soluprot] 3/3 downloading 64-bit USEARCH v11 (public domain)"
mkdir -p "$INSTALL_DIR/usearch"
curl -fsSL https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/main/bin/usearch11.0.667_i86linux64 \
    -o "$INSTALL_DIR/usearch/usearch11.0.667_i86linux64"
chmod +x "$INSTALL_DIR/usearch/usearch11.0.667_i86linux64"

# Point soluprot.py's built-in tool paths at the absolute locations above.
sed -i "s#^    _USEARCH = .*#    _USEARCH = '$INSTALL_DIR/usearch/usearch11.0.667_i86linux64'#" "$SOLUPROT_PATH/soluprot.py"
sed -i "s#^    _TMHMM = .*#    _TMHMM = '$INSTALL_DIR/tmhmm/bin/tmhmm'#" "$SOLUPROT_PATH/soluprot.py"

cat <<EOF

==============================================================================
SoluProt env + code + USEARCH installed.
  env:          $SOLUPROT_ENV
  SoluProt:     $SOLUPROT_PATH
  USEARCH:      $INSTALL_DIR/usearch/usearch11.0.667_i86linux64

STILL MANUAL — TMHMM 2.0 (academic license):
  Obtain from https://services.healthtech.dtu.dk/ (TMHMM 2.0) or
  https://git.loschmidt.cz/misc/tmhmm . Unpack so that
  $INSTALL_DIR/tmhmm/bin/tmhmm exists, then fix the perl shebang:
    sed -i '1 s|^#!.*perl.*|#!/usr/bin/env perl|' \\
      "$INSTALL_DIR"/tmhmm/bin/tmhmm "$INSTALL_DIR"/tmhmm/bin/tmhmmformat.pl
  (Its 64-bit binary is bin/decodeanhmm.Linux_x86_64.)
  Without TMHMM, SoluProt can still run via --no_tmhmm + the notmhmm model,
  but run_soluprot.sh does not pass --no_tmhmm, so install TMHMM for default use.

Finally, set in tps_eval/paths.sh:
  SOLUPROT_ENV="$SOLUPROT_ENV"
  SOLUPROT_PATH="$SOLUPROT_PATH"

Self-test:
  conda activate $SOLUPROT_ENV
  cd "$SOLUPROT_PATH"
  python soluprot.py --i_fa data/test.fa --o_csv /tmp/test.csv --tmp_dir /tmp
  diff /tmp/test.csv data/test.csv   # should be identical
==============================================================================
EOF
