#!/bin/bash
# Build the interactive MARTS-DB natural-bands HTML dashboard. Standalone reporting
# utility (NOT part of the eval / submit-all pipeline); runs locally on a login node
# or laptop. Pure Python stdlib — no conda env required. See src/dashboard/ for logic.

USAGE="[--designs <merged.csv|metrics_dir>] [--demo] [--bands <json...>] [--output <out.html>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Renders the committed MARTS-DB reference bands (src/reference_stats/*.json) as a"
    echo "self-contained interactive HTML page: per-metric percentile envelopes, stratifiable"
    echo "by substrate / first-cyclization / domain architecture, across the ESMFold, AF3 and"
    echo "Boltz-2 (holo) structure sources. Optionally overlays a generated-design batch."
    echo
    echo "Arguments:"
    echo "  --designs PATH  A design batch to overlay: a merged CSV or a directory of the"
    echo "                  pipeline's *_<tool>.csv outputs (matched to bands by column name)."
    echo "  --demo          Overlay a synthetic demo batch (ignored if --designs is given)."
    echo "  --bands JSON…   Band JSON paths (default: the committed esmfold/af3/boltz2 JSONs)."
    echo "  --output PATH   Output HTML (default: data/dashboard/marts_db_dashboard.html)."
    echo "  -h, --help      Show this help and exit."
}

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    Help
    exit 0
fi

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.." || exit 1

PY=${PYTHON:-}
if [[ -z "$PY" ]]; then
    if command -v python3 >/dev/null 2>&1; then PY=python3; else PY=python; fi
fi
"$PY" src/dashboard/build_dashboard.py "$@"
