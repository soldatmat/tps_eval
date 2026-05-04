#!/bin/bash
# Karolina-specific environment setup for tps_eval.
#
# Functionally equivalent to the project-root setup.sh (creates the conda env
# from the same channels and packages, then `pip install` torch + fair-esm and
# `pip install -e .` for this repo). Differences are purely Karolina adaptations:
#
#   * Loads the cluster's Anaconda3 module instead of assuming conda is on PATH.
#   * Creates the env at a path-based location on project storage
#     (/mnt/proj2/fta-26-15/.conda/envs/tps_eval) rather than under $HOME — the
#     user's $HOME quota (~23 GB) cannot fit a full env plus repodata cache.
#   * Redirects CONDA_PKGS_DIRS, PIP_CACHE_DIR, and XDG_CACHE_HOME to project
#     storage for the same quota reason.
#   * Removes any stale partial env in $HOME and any existing target env so the
#     script is idempotent — re-running always yields a clean install.
#
# Runnable from any directory; all paths inside are absolute.
set -euo pipefail

PROJECT_ENV=/mnt/proj2/fta-26-15/.conda/envs/tps_eval
PARTIAL_ENV=/home/soldatmat/.conda/envs/tps_eval
TPS_EVAL_REPO=/mnt/proj2/fta-26-15/documents/tps_eval

# Home quota on Karolina is ~23 GB and easily fills with conda repodata + pkgs.
# Redirect every cache that conda/pip might use to project storage.
export CONDA_PKGS_DIRS=/mnt/proj2/fta-26-15/.conda/pkgs
export PIP_CACHE_DIR=/mnt/proj2/fta-26-15/.cache/pip
export XDG_CACHE_HOME=/mnt/proj2/fta-26-15/.cache
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME"

echo "[$(date '+%H:%M:%S')] Loading Anaconda3 module"
module purge >/dev/null 2>&1 || true
ml Anaconda3

eval "$(conda shell.bash hook)"

echo "[$(date '+%H:%M:%S')] Deleting partial env at $PARTIAL_ENV"
rm -rf "$PARTIAL_ENV"

if [[ -d "$PROJECT_ENV" ]]; then
    echo "[$(date '+%H:%M:%S')] Existing env at $PROJECT_ENV — removing for clean install"
    rm -rf "$PROJECT_ENV"
fi

echo "[$(date '+%H:%M:%S')] Creating env at $PROJECT_ENV (this may take a while)"
conda create -p "$PROJECT_ENV" \
    -c conda-forge -c bioconda -c schrodinger \
    python biopython pandas requests tqdm openbabel foldseek pymol-bundle \
    -y

echo "[$(date '+%H:%M:%S')] Activating env"
conda activate "$PROJECT_ENV"

echo "[$(date '+%H:%M:%S')] Active python: $(which python)"
echo "[$(date '+%H:%M:%S')] pip install torch"
pip install torch

echo "[$(date '+%H:%M:%S')] pip install fair-esm"
pip install fair-esm

echo "[$(date '+%H:%M:%S')] pip install -e $TPS_EVAL_REPO"
cd "$TPS_EVAL_REPO"
pip install -e .

echo "[$(date '+%H:%M:%S')] DONE — env installed at $PROJECT_ENV"
