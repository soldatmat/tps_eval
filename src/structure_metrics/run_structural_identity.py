from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from foldseek.structure_alignment import main as _structure_alignment  # noqa: E402

# foldseek best-hit columns -> pipeline metric names (keyed by ID), the structural
# analog of max_sequence_identity: per generated structure, the best (max) TM-score
# and lddt to the nearest known-TPS structure, plus which reference it matched.
_RENAME = {
    "query": "ID",
    "max_alntmscore": "structural_tmscore_to_known",
    "max_alntmscore_target": "structural_tmscore_to_known_hit",
    "max_lddt": "structural_lddt_to_known",
    "max_lddt_target": "structural_lddt_to_known_hit",
}


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_structural_identity.csv")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Foldseek structural similarity of each generated structure to the "
        "nearest known-TPS reference structure (max TM-score / lddt). The structural "
        "analog of max_sequence_identity. Writes a CSV keyed by ID."
    )
    parser.add_argument("structs_dir", help="Directory of generated structures (.pdb/.cif).")
    parser.add_argument("known_structs_dir", help="Directory of known-TPS reference structures.")
    parser.add_argument("--save_path", default=None,
                        help="Output CSV path (default: <structs_dir>_structural_identity.csv).")
    args = parser.parse_args()

    tmp = tempfile.mkdtemp(prefix="structural_identity_")
    try:
        _structure_alignment(argparse.Namespace(
            structures_root=args.structs_dir,
            known_structures_root=args.known_structs_dir,
            output_root=tmp,
            store_intermediate_results=False,
            random_run_id=False,
        ))
        scores = pd.read_csv(os.path.join(tmp, "structure_alignment_scores.csv"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    cols = [c for c in _RENAME if c in scores.columns]
    scores = scores[cols].rename(columns=_RENAME)
    # Strip any text after whitespace in IDs, matching the other metrics' join key.
    if "ID" in scores.columns:
        scores["ID"] = scores["ID"].astype(str).map(lambda x: x.split(" ", 1)[0])

    save_path = args.save_path or _default_save_path(args.structs_dir)
    scores.to_csv(save_path, index=False)
    print(f"Wrote {len(scores)} rows to {save_path}")


if __name__ == "__main__":
    main()
