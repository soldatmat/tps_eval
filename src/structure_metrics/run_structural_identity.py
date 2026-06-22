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


def _default_topk_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(
        os.path.dirname(d), os.path.basename(d) + "_structural_identity_topk.csv"
    )


def _stem(name: str) -> str:
    """Target structure stem: drop any path and a trailing .pdb/.cif/.gz extension."""
    base = os.path.basename(str(name))
    for ext in (".pdb.gz", ".cif.gz", ".pdb", ".cif", ".ent"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def _write_topk(raw_csv_path: str, top_k: int, save_path: str) -> None:
    """Top-k targets per query by alntmscore (LARGER = closer), from foldseek hits.

    Tidy CSV with columns query_id,rank,neighbour_id,score where score is the
    foldseek alntmscore (TM-score) and neighbour_id is the target structure stem.
    foldseek returns many hits per query; we keep the k highest by alntmscore.
    Self-hits (target stem == query stem) are excluded so the same contract holds
    when a structure set is searched against itself (leave-one-out).
    """
    hits = pd.read_csv(raw_csv_path)
    rows = []
    for query, grp in hits.groupby("query", sort=True):
        query_id = str(query).split(" ", 1)[0]
        query_stem = _stem(query)
        grp = grp.copy()
        grp["__nid"] = grp["target"].map(_stem)
        # Exclude self-hits (handles self-search leave-one-out).
        grp = grp[grp["__nid"] != query_stem]
        # Highest alntmscore first; stable sort keeps foldseek's order on ties.
        grp = grp.sort_values("alntmscore", ascending=False, kind="stable")
        for rank, (_, row) in enumerate(grp.head(top_k).iterrows(), start=1):
            rows.append(
                {
                    "query_id": query_id,
                    "rank": rank,
                    "neighbour_id": row["__nid"],
                    "score": float(row["alntmscore"]),
                }
            )
    pd.DataFrame(rows, columns=["query_id", "rank", "neighbour_id", "score"]).to_csv(
        save_path, index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Foldseek structural similarity of each generated structure to the "
        "nearest known-TPS reference structure (max TM-score / lddt). The structural "
        "analog of max_sequence_identity. Writes a CSV keyed by ID. With --top_k N, "
        "additionally writes <structs_dir>_structural_identity_topk.csv (columns "
        "query_id,rank,neighbour_id,score) where score is the foldseek alntmscore "
        "(TM-score; LARGER = closer) and neighbour_id is the target structure stem. "
        "Default single-best output is unchanged."
    )
    parser.add_argument("structs_dir", help="Directory of generated structures (.pdb/.cif).")
    parser.add_argument("known_structs_dir", help="Directory of known-TPS reference structures.")
    parser.add_argument("--save_path", default=None,
                        help="Output CSV path (default: <structs_dir>_structural_identity.csv).")
    parser.add_argument("--top_k", type=int, default=None,
                        help="If >= 1, also emit the top-k nearest reference structures per query.")
    parser.add_argument("--self_mode", action="store_true", default=False,
                        help="Searching a structure set against itself: drop self-hits "
                        "(target stem == query stem) before the single-best reduction, so "
                        "each query's best hit is its nearest OTHER neighbour (leave-one-out) "
                        "instead of the trivial self-match TM~1.0. The top-k path already "
                        "excludes self-hits unconditionally.")
    args = parser.parse_args()

    want_topk = args.top_k is not None and args.top_k >= 1

    tmp = tempfile.mkdtemp(prefix="structural_identity_")
    try:
        _structure_alignment(argparse.Namespace(
            structures_root=args.structs_dir,
            known_structures_root=args.known_structs_dir,
            output_root=tmp,
            # Keep the raw per-hit table so top-k can be derived from it.
            store_intermediate_results=want_topk,
            random_run_id=False,
            exclude_self=args.self_mode,
        ))
        scores = pd.read_csv(os.path.join(tmp, "structure_alignment_scores.csv"))

        if want_topk:
            raw_csv = os.path.join(tmp, "structure_alignments.csv")
            topk_save_path = _default_topk_save_path(args.structs_dir)
            _write_topk(raw_csv, args.top_k, topk_save_path)
            print(f"Wrote top-{args.top_k} neighbours to {topk_save_path}")
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
