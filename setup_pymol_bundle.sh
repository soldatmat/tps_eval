#!/bin/bash
# Optional add-on to setup.sh: replace pymol-open-source with the licensed
# Schrodinger pymol-bundle (Incentive PyMOL). Run this AFTER the default
# tps_eval env is already installed via setup.sh.
#
# Why one might want this: Incentive PyMOL has a higher-quality ray-tracer
# (better shadows/AO/transparency) and bundled plugins like APBS. For the
# scripted plot_domains / plot_residue_similarity workflow the open-source
# build is functionally identical — this swap is purely a render-quality /
# plugin choice.
#
# Without a valid license file, pymol-bundle renders carry a Schrodinger
# watermark. PyMOL Incentive looks for the license file in (in order):
#   1. $PYMOL_LICENSE_FILE (env var)
#   2. ~/.pymol/license.lic
#   3. <conda_prefix>/share/pymol/license.lic
#   4. /etc/pymol/license.lic
# Pass --license_path to copy a license into <conda_prefix>/share/pymol/
# so it's auto-discovered without touching $HOME or env vars.

USAGE="[--license_path <path>]"

Help() {
    echo "Usage: $0 $USAGE"
    echo
    echo "Replaces pymol-open-source in the existing tps_eval conda env with"
    echo "Schrodinger pymol-bundle (Incentive PyMOL)."
    echo
    echo "Arguments:"
    echo "  --license_path  Path to a PyMOL Incentive .lic file. If provided,"
    echo "                  copied into <conda_prefix>/share/pymol/license.lic"
    echo "                  so PyMOL discovers it automatically."
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

eval "$(conda shell.bash hook)"
conda activate tps_eval

echo "[$(date '+%H:%M:%S')] Removing pymol-open-source"
conda remove -n tps_eval -y pymol-open-source || true

echo "[$(date '+%H:%M:%S')] Installing pymol-bundle (with catch2=3.13 ABI pin)"
conda install -n tps_eval -c conda-forge -c schrodinger -y \
    pymol-bundle "catch2=3.13"

if [[ -n "$license_path" ]]; then
    if [[ ! -f "$license_path" ]]; then
        echo "ERROR: license file not found at $license_path" >&2
        exit 1
    fi
    dest="$CONDA_PREFIX/share/pymol/license.lic"
    mkdir -p "$(dirname "$dest")"
    cp "$license_path" "$dest"
    echo "[$(date '+%H:%M:%S')] Installed license at $dest"
else
    echo "[$(date '+%H:%M:%S')] No --license_path given. To activate the license later, drop your .lic file at one of:"
    echo "  - \$PYMOL_LICENSE_FILE (set as env var)"
    echo "  - ~/.pymol/license.lic"
    echo "  - $CONDA_PREFIX/share/pymol/license.lic"
fi

echo "[$(date '+%H:%M:%S')] DONE"
