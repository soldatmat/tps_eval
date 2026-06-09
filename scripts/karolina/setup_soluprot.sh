#!/bin/bash
# Karolina-specific SoluProt setup. Same end result as scripts/setup_soluprot.sh
# (py3.7 env from soluprot_environment.yml + SoluProt standalone + 64-bit USEARCH;
# TMHMM manual), with Karolina adaptations:
#   * loads the Anaconda3 module
#   * installs under project storage (/mnt/proj2/fta-26-15), never $HOME (~23 GB quota)
#   * redirects conda/pip/XDG caches to project storage
#   * creates the env at a path and appends that dir to conda envs_dirs so
#     `conda activate soluprot` BY NAME works (run_soluprot.sh activates by name)
#
# Runnable from any directory; all paths absolute.
set -euo pipefail

PROJECT=/mnt/proj2/fta-26-15
INSTALL_DIR="$PROJECT/documents"            # soluprot/, usearch/, tmhmm/ go here
SOLUPROT_ENV_DIR="$PROJECT/.conda/envs/soluprot"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # tps_eval/scripts/karolina

export CONDA_PKGS_DIRS="$PROJECT/.conda/pkgs"
export PIP_CACHE_DIR="$PROJECT/.cache/pip"
export XDG_CACHE_HOME="$PROJECT/.cache"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME" "$INSTALL_DIR"

echo "[$(date '+%H:%M:%S')] Loading Anaconda3"
module purge >/dev/null 2>&1 || true
ml Anaconda3
eval "$(conda shell.bash hook)"

echo "[$(date '+%H:%M:%S')] 1/3 creating soluprot env at $SOLUPROT_ENV_DIR"
[ -d "$SOLUPROT_ENV_DIR" ] && rm -rf "$SOLUPROT_ENV_DIR"
conda env create -p "$SOLUPROT_ENV_DIR" -f "$SCRIPT_DIR/../soluprot_environment.yml"
# So `conda activate soluprot` by name resolves (run_soluprot.sh uses the name):
conda config --append envs_dirs "$PROJECT/.conda/envs" 2>/dev/null || true

echo "[$(date '+%H:%M:%S')] 2/3 downloading SoluProt standalone"
cd "$INSTALL_DIR"
curl -fsSL 'https://loschmidt.chemi.muni.cz/soluprot/?f=soluprot.zip' -o soluprot.zip
unzip -qo soluprot.zip
rm -f soluprot.zip
SOLUPROT_PATH="$(dirname "$(find "$INSTALL_DIR" -maxdepth 3 -name soluprot.py | head -1)")"

echo "[$(date '+%H:%M:%S')] 3/3 downloading 64-bit USEARCH v11"
mkdir -p "$INSTALL_DIR/usearch"
curl -fsSL https://raw.githubusercontent.com/rcedgar/usearch_old_binaries/main/bin/usearch11.0.667_i86linux64 \
    -o "$INSTALL_DIR/usearch/usearch11.0.667_i86linux64"
chmod +x "$INSTALL_DIR/usearch/usearch11.0.667_i86linux64"

sed -i "s#^    _USEARCH = .*#    _USEARCH = '$INSTALL_DIR/usearch/usearch11.0.667_i86linux64'#" "$SOLUPROT_PATH/soluprot.py"
sed -i "s#^    _TMHMM = .*#    _TMHMM = '$INSTALL_DIR/tmhmm/bin/tmhmm'#" "$SOLUPROT_PATH/soluprot.py"

cat <<EOF

[$(date '+%H:%M:%S')] DONE (env + code + USEARCH).
  env:      $SOLUPROT_ENV_DIR  (activate by name: 'soluprot')
  SoluProt: $SOLUPROT_PATH
  USEARCH:  $INSTALL_DIR/usearch/usearch11.0.667_i86linux64

TMHMM 2.0 is still MANUAL (academic license) — unpack to $INSTALL_DIR/tmhmm so
$INSTALL_DIR/tmhmm/bin/tmhmm exists, use bin/decodeanhmm.Linux_x86_64, and fix
the perl shebang in bin/tmhmm + bin/tmhmmformat.pl (#!/usr/bin/env perl).

Then set in $PROJECT/documents/tps_eval/paths.sh:
  SOLUPROT_ENV="soluprot"
  SOLUPROT_PATH="$SOLUPROT_PATH"

NOTE: tps_eval-env tools on Karolina need the libstdc++ LD_LIBRARY_PATH fix; the
soluprot env (py3.7) does not hit it, but run on a compute node if in doubt.
EOF
