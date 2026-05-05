#!/bin/bash
# Karolina-specific add-on to scripts/karolina/setup.sh: replace
# pymol-open-source with the licensed Schrodinger pymol-bundle (Incentive
# PyMOL). Run AFTER scripts/karolina/setup.sh has installed the env at
# $PROJECT_ENV.
#
# See ../../setup_pymol_bundle.sh for the rationale. This Karolina variant:
#   * Loads the Anaconda3 module instead of assuming conda is on PATH.
#   * Targets the path-based env at /mnt/proj2/fta-26-15/.conda/envs/tps_eval.
#   * Redirects CONDA_PKGS_DIRS / PIP_CACHE_DIR / XDG_CACHE_HOME to project
#     storage so the user's ~23 GB $HOME quota is not consumed.
#
# Without a valid license file, pymol-bundle renders carry a Schrodinger
# watermark. PyMOL Incentive looks for the license file in (in order):
#   1. $PYMOL_LICENSE_FILE (env var)
#   2. ~/.pymol/license.lic
#   3. <conda_prefix>/share/pymol/license.lic
#   4. /etc/pymol/license.lic
# Pass --license_path to copy a license into <conda_prefix>/share/pymol/.

USAGE="[--license_path <path>]"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --license_path  Path to a PyMOL Incentive .lic file (optional)"
    echo "  -h, --help      Show this help message and exit"
    echo
}

license_path=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --license_path) license_path="$2"; shift 2 ;;
        -h|--help) Help; exit 0 ;;
        *) echo "Unknown option: $1"; Help; exit 1 ;;
    esac
done

set -euo pipefail

PROJECT_ENV=/mnt/proj2/fta-26-15/.conda/envs/tps_eval

export CONDA_PKGS_DIRS=/mnt/proj2/fta-26-15/.conda/pkgs
export PIP_CACHE_DIR=/mnt/proj2/fta-26-15/.cache/pip
export XDG_CACHE_HOME=/mnt/proj2/fta-26-15/.cache
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME"

echo "[$(date '+%H:%M:%S')] Loading Anaconda3 module"
module purge >/dev/null 2>&1 || true
ml Anaconda3

eval "$(conda shell.bash hook)"

if [[ ! -d "$PROJECT_ENV" ]]; then
    echo "ERROR: env not found at $PROJECT_ENV. Run scripts/karolina/setup.sh first." >&2
    exit 1
fi

conda activate "$PROJECT_ENV"

echo "[$(date '+%H:%M:%S')] Removing pymol-open-source"
conda remove -p "$PROJECT_ENV" -y pymol-open-source || true

echo "[$(date '+%H:%M:%S')] Installing pymol-bundle (with catch2=3.13 ABI pin)"
conda install -p "$PROJECT_ENV" -c conda-forge -c schrodinger -y \
    pymol-bundle "catch2=3.13"

if [[ -n "$license_path" ]]; then
    if [[ ! -f "$license_path" ]]; then
        echo "ERROR: license file not found at $license_path" >&2
        exit 1
    fi
    dest="$PROJECT_ENV/share/pymol/license.lic"
    mkdir -p "$(dirname "$dest")"
    cp "$license_path" "$dest"
    echo "[$(date '+%H:%M:%S')] Installed license at $dest"
else
    echo "[$(date '+%H:%M:%S')] No --license_path given. To activate the license later, drop your .lic file at one of:"
    echo "  - \$PYMOL_LICENSE_FILE (set as env var)"
    echo "  - ~/.pymol/license.lic"
    echo "  - $PROJECT_ENV/share/pymol/license.lic"
fi

echo "[$(date '+%H:%M:%S')] DONE"
