#!/bin/bash
# Setup for CataPro (enzyme-kinetics predictor), used by run_catapro.sh. CataPro is
# vendored as a git submodule at vendor/CataPro (its trained 10-fold heads —
# kcat_models/ Km_models/ act_models/ — ship inside the submodule). This script
# automates the two parts that are NOT in the repo:
#   1. the conda env from scripts/setup/catapro_environment.yml
#   2. the two HuggingFace backbones CataPro needs at inference time, placed inside
#      vendor/CataPro/models/ where predict.py expects them:
#        - Rostlab/prot_t5_xl_uniref50      (enzyme features; ~large, several GB)
#        - laituan245/molt5-base-smiles2caption (substrate features)
# These backbones are NOT committed (large; per the repo's DB/weights convention).
#
# Usage:  ./setup_catapro.sh
# Assumes `conda` is already available and the CataPro submodule is checked out
# (git submodule update --init vendor/CataPro).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"   # scripts/setup -> repo root
CATAPRO_ENV="${CATAPRO_ENV:-catapro}"
MODELS_DIR="$REPO_ROOT/vendor/CataPro/models"

if [ ! -f "$REPO_ROOT/vendor/CataPro/inference/predict.py" ]; then
    echo "ERROR: vendor/CataPro not checked out. Run: git submodule update --init vendor/CataPro" >&2
    exit 1
fi

echo "[setup_catapro] 1/2 creating conda env '$CATAPRO_ENV'"
conda env create -n "$CATAPRO_ENV" -f "$SCRIPT_DIR/catapro_environment.yml"

echo "[setup_catapro] 2/2 downloading ProtT5 + MolT5 backbones into $MODELS_DIR"
mkdir -p "$MODELS_DIR"
eval "$(conda shell.bash hook)"
conda activate "$CATAPRO_ENV"
python - "$MODELS_DIR" <<'PY'
import sys
from huggingface_hub import snapshot_download

models_dir = sys.argv[1]
# Grab config + tokenizer + PyTorch weights only (skip TF/Flax/rust duplicates).
allow = ["*.json", "*.bin", "*.model", "*.txt", "*.safetensors"]
for repo, subdir in [
    ("Rostlab/prot_t5_xl_uniref50", "prot_t5_xl_uniref50"),
    ("laituan245/molt5-base-smiles2caption", "molt5-base-smiles2caption"),
]:
    dest = f"{models_dir}/{subdir}"
    print(f"  downloading {repo} -> {dest}")
    snapshot_download(repo_id=repo, local_dir=dest,
                      local_dir_use_symlinks=False, allow_patterns=allow)
print("  done.")
PY

cat <<EOF

==============================================================================
CataPro env + backbones installed.
  env:      $CATAPRO_ENV
  models:   $MODELS_DIR (prot_t5_xl_uniref50, molt5-base-smiles2caption,
            plus the in-repo kcat_models/ Km_models/ act_models/)

Set in tps_eval/paths.sh:
  CATAPRO_ENV="$CATAPRO_ENV"

Self-test (on a GPU node):
  conda activate $CATAPRO_ENV
  cd "$REPO_ROOT"
  sh scripts/tool_wrappers/run_catapro.sh --fasta_path <some.fasta> --target_substrate FPP --device cuda
  # writes <some>_catapro.csv keyed by ID (catapro_kcat, catapro_km, catapro_kcat_km)
==============================================================================
EOF
