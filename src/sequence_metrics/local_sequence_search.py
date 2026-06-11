from __future__ import annotations

"""Fast LOCAL (BLAST-style) sequence identity/similarity search + top-k neighbours.

This is the LOCAL counterpart of ``max_sequence_identity`` (which stays as the GLOBAL,
full-length Biopython novelty metric). Instead of Needleman-Wunsch over every query x
reference pair, it delegates to a fast, backend-pluggable local aligner:

  * **mmseqs2** (DEFAULT) — ``mmseqs easy-search`` (Smith-Waterman-ish local hits).
  * **diamond**           — ``diamond blastp`` (local hits).

For each query it reports the BEST-HIT local identity (+ similarity + coverage), and,
when ``top_k`` is given, the k closest reference neighbours as a tidy CSV with EXACTLY
``query_id,rank,neighbour_id,score`` (score = identity %, LARGER = closer) — the same
contract as the existing ``_topk.csv`` files, so the k-NN / SDR sequence space can
consume it directly.

Why a separate tool (not a replacement): the global metric answers "how novel is this
design overall" (full-length identity), which matters for the novelty story; the local
metric answers "what's the closest local relative and how close" in SECONDS rather than
hours, and is what the k-NN calibration actually needs. The Biopython self all-vs-all on
the 1195-seq MARTS reference took >4h and TIMED OUT; the local backends do it in
seconds-to-minutes.

Per-backend output-field mapping (mapped CONSISTENTLY so callers don't care which backend
produced the CSV):

    metric column                  mmseqs2 field        diamond field
    -----------------------------  -------------------  ---------------
    local_sequence_identity        fident (* 100)       pident
    local_sequence_similarity      NaN (see note)       ppos
    local_coverage                 qcov   (* 100)       qcovhsp

Note on similarity: DIAMOND emits ``ppos`` (% positive-scoring columns) directly. MMseqs2
``easy-search`` has NO positives/similarity field in its ``--format-output`` vocabulary
(it exposes fident/pident/qcov/tcov/alnlen/mismatch/gapopen/... but nothing equivalent to
BLAST ``ppos``). Rather than fake it, the mmseqs2 backend emits ``local_sequence_similarity``
as NaN. Identity is the PRIMARY field and is available from BOTH backends; similarity is a
DIAMOND-only extra. Output CSV is keyed by ``ID``: ``<input>_local_sequence_search.csv``.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers  # noqa: E402


BACKENDS = ("mmseqs2", "diamond")

METRIC_COLUMNS = [
    "ID",
    "local_sequence_identity",
    "local_sequence_identity_hit",
    "local_sequence_similarity",
    "local_sequence_similarity_hit",
    "local_coverage",
]

TOPK_COLUMNS = ["query_id", "rank", "neighbour_id", "score"]

# Internal normalized hit columns produced by each backend's parser.
# identity/similarity/coverage are PERCENT in [0, 100]; score = identity %.
_HIT_COLUMNS = ["qseqid", "sseqid", "identity", "similarity", "coverage", "score"]


def _first_token(value: str) -> str:
    return str(value).split(" ", 1)[0].split("\t", 1)[0]


# ---------------------------------------------------------------------------
# Backend: MMseqs2 (easy-search)
# ---------------------------------------------------------------------------
# easy-search --format-output fields. fident = fraction identity in [0,1];
# qcov = query coverage in [0,1]. There is no positives/similarity field, so
# similarity is filled with NaN downstream for this backend.
_MMSEQS_FORMAT = "query,target,fident,qcov,bits,evalue"


def _run_mmseqs(fasta_path: str, ref_fasta: str, work_dir: str, *,
                top_k: int, threads: int, sensitivity: Optional[str]) -> pd.DataFrame:
    out_tsv = os.path.join(work_dir, "hits.m8")
    tmp_dir = os.path.join(work_dir, "mmseqs_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = [
        "mmseqs", "easy-search",
        fasta_path, ref_fasta, out_tsv, tmp_dir,
        "--format-output", _MMSEQS_FORMAT,
        # Keep enough hits per query to populate top-k neighbours.
        "--max-seqs", str(max(top_k, 1) * 5 + 50),
        "--threads", str(threads),
    ]
    if sensitivity is not None:
        cmd += ["-s", str(sensitivity)]
    _run(cmd)

    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        return pd.DataFrame(columns=_HIT_COLUMNS)

    raw = pd.read_csv(
        out_tsv, sep="\t", header=None,
        names=["query", "target", "fident", "qcov", "bits", "evalue"],
    )
    hits = pd.DataFrame({
        "qseqid": raw["query"].astype(str).map(_first_token),
        "sseqid": raw["target"].astype(str).map(_first_token),
        "identity": raw["fident"].astype(float) * 100.0,
        # MMseqs2 has no positives field -> similarity unavailable.
        "similarity": np.nan,
        "coverage": raw["qcov"].astype(float) * 100.0,
        # Rank/best-hit selection uses bitscore; score column = identity %.
        "_bits": raw["bits"].astype(float),
    })
    hits["score"] = hits["identity"]
    return hits


# ---------------------------------------------------------------------------
# Backend: DIAMOND (blastp)
# ---------------------------------------------------------------------------
# DIAMOND outfmt 6 fields. pident/ppos/qcovhsp are already PERCENT in [0,100].
_DIAMOND_OUTFMT = ["qseqid", "sseqid", "pident", "ppos", "qcovhsp", "bitscore"]


def _run_diamond(fasta_path: str, ref_fasta: str, work_dir: str, *,
                 top_k: int, threads: int, sensitivity: Optional[str]) -> pd.DataFrame:
    db = os.path.join(work_dir, "ref_db")
    _run(["diamond", "makedb", "--in", ref_fasta, "--db", db,
          "--threads", str(threads), "--quiet"])

    out_tsv = os.path.join(work_dir, "hits.tsv")
    cmd = [
        "diamond", "blastp",
        "--db", db,
        "--query", fasta_path,
        "--out", out_tsv,
        "--outfmt", "6", *_DIAMOND_OUTFMT,
        "--max-target-seqs", str(max(top_k, 1) * 5 + 50),
        "--threads", str(threads),
        "--quiet",
    ]
    cmd += [f"--{sensitivity}"] if sensitivity else ["--very-sensitive"]
    _run(cmd)

    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        return pd.DataFrame(columns=_HIT_COLUMNS)

    raw = pd.read_csv(out_tsv, sep="\t", header=None, names=_DIAMOND_OUTFMT)
    hits = pd.DataFrame({
        "qseqid": raw["qseqid"].astype(str).map(_first_token),
        "sseqid": raw["sseqid"].astype(str).map(_first_token),
        "identity": raw["pident"].astype(float),
        "similarity": raw["ppos"].astype(float),
        "coverage": raw["qcovhsp"].astype(float),
        "_bits": raw["bitscore"].astype(float),
    })
    hits["score"] = hits["identity"]
    return hits


_BACKEND_RUNNERS = {
    "mmseqs2": _run_mmseqs,
    "diamond": _run_diamond,
}


def _run(cmd: List[str]) -> None:
    print("Running:", " ".join(str(c) for c in cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end="")
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)


# ---------------------------------------------------------------------------
# Reductions
# ---------------------------------------------------------------------------
def _best_hits(hits: pd.DataFrame, *, self_mode: bool) -> pd.DataFrame:
    """Best (max-bitscore) hit per query. In self mode the query's own id is excluded."""
    if hits.empty:
        return hits
    if self_mode:
        hits = hits[hits["qseqid"] != hits["sseqid"]]
    if hits.empty:
        return hits
    # Highest bitscore first; deterministic tie-break on identity then subject id.
    hits = hits.sort_values(
        ["_bits", "identity", "sseqid"], ascending=[False, False, True]
    )
    return hits.groupby("qseqid", as_index=False, sort=False).first()


def _topk_neighbours(hits: pd.DataFrame, query_ids: List[str], top_k: int, *,
                     self_mode: bool) -> pd.DataFrame:
    """Tidy top-k neighbours per query, ranked by score (identity %, descending)."""
    rows: List[dict] = []
    if not hits.empty:
        if self_mode:
            hits = hits[hits["qseqid"] != hits["sseqid"]]
        # Deduplicate to one row per (query, neighbour) keeping the best bitscore,
        # then rank by score (identity %), tie-break on bitscore then subject id.
        hits = hits.sort_values(
            ["qseqid", "_bits"], ascending=[True, False]
        ).drop_duplicates(["qseqid", "sseqid"], keep="first")
        by_query = {q: g for q, g in hits.groupby("qseqid", sort=False)}
    else:
        by_query = {}

    for qid in query_ids:
        group = by_query.get(qid)
        if group is None:
            continue
        group = group.sort_values(
            ["score", "_bits", "sseqid"], ascending=[False, False, True]
        ).head(top_k)
        for rank, (_, row) in enumerate(group.iterrows(), start=1):
            rows.append({
                "query_id": qid,
                "rank": rank,
                "neighbour_id": row["sseqid"],
                "score": float(row["score"]),
            })
    return pd.DataFrame(rows, columns=TOPK_COLUMNS)


# ---------------------------------------------------------------------------
# Public entry
# ---------------------------------------------------------------------------
def _default_save_path(fasta_path: str) -> str:
    base, _ext = os.path.splitext(fasta_path)
    return f"{base}_local_sequence_search.csv"


def _default_topk_save_path(fasta_path: str) -> str:
    base, _ext = os.path.splitext(fasta_path)
    return f"{base}_local_sequence_search_topk.csv"


def local_sequence_search(
    fasta_path: str,
    *,
    train_path: Optional[str] = None,
    backend: str = "mmseqs2",
    top_k: Optional[int] = None,
    threads: int = 4,
    sensitivity: Optional[str] = None,
    save_path: Optional[str] = None,
    topk_save_path: Optional[str] = None,
) -> pd.DataFrame:
    """Run a local identity/similarity search of ``fasta_path`` vs a reference.

    Args:
        fasta_path: Query FASTA.
        train_path: Reference FASTA. If omitted, SELF mode (each query's best hit /
            neighbours exclude itself) — mirrors max_sequence_identity's self mode.
        backend: ``mmseqs2`` (default) or ``diamond``.
        top_k: If >= 1, also write the tidy top-k neighbours CSV.
        threads: Worker threads passed to the backend.
        sensitivity: Backend sensitivity knob (mmseqs2 ``-s`` value, e.g. 7.5; diamond
            flag name, e.g. very-sensitive). None -> backend default.
        save_path / topk_save_path: Output paths (defaults derived from fasta_path).
    """
    if backend not in BACKENDS:
        raise ValueError(f"Unknown backend {backend!r}; choose from {BACKENDS}.")

    self_mode = train_path is None
    ref_fasta = fasta_path if self_mode else train_path

    query_identifiers, _seqs = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    query_ids = [_first_token(i) for i in query_identifiers]

    runner = _BACKEND_RUNNERS[backend]
    effective_top_k = top_k if (top_k and top_k >= 1) else 1

    work_dir = tempfile.mkdtemp(prefix=f"local_seq_search_{backend}_")
    try:
        hits = runner(
            fasta_path, ref_fasta, work_dir,
            top_k=effective_top_k, threads=threads, sensitivity=sensitivity,
        )
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    best = _best_hits(hits, self_mode=self_mode)
    best_by_query = {row["qseqid"]: row for _, row in best.iterrows()} if not best.empty else {}

    rows: List[dict] = []
    for qid in query_ids:
        b = best_by_query.get(qid)
        if b is None:
            rows.append({
                "ID": qid,
                "local_sequence_identity": np.nan,
                "local_sequence_identity_hit": "",
                "local_sequence_similarity": np.nan,
                "local_sequence_similarity_hit": "",
                "local_coverage": np.nan,
            })
        else:
            rows.append({
                "ID": qid,
                "local_sequence_identity": float(b["identity"]),
                "local_sequence_identity_hit": b["sseqid"],
                "local_sequence_similarity": float(b["similarity"]),
                # similarity hit == identity hit (single best hit per query).
                "local_sequence_similarity_hit": (
                    b["sseqid"] if not pd.isna(b["similarity"]) else ""
                ),
                "local_coverage": float(b["coverage"]),
            })

    df = pd.DataFrame(rows)[METRIC_COLUMNS]
    if save_path is None:
        save_path = _default_save_path(fasta_path)
    df.to_csv(save_path, index=False)
    n_with_hits = int(df["local_sequence_identity"].notna().sum())
    print(f"[{backend}] queries: {len(df)}  with >=1 hit: {n_with_hits}  "
          f"-> {save_path}")

    if top_k and top_k >= 1:
        topk_df = _topk_neighbours(hits, query_ids, top_k, self_mode=self_mode)
        if topk_save_path is None:
            topk_save_path = _default_topk_save_path(fasta_path)
        topk_df.to_csv(topk_save_path, index=False)
        print(f"[{backend}] top-{top_k} neighbours: {len(topk_df)} rows -> {topk_save_path}")

    return df
