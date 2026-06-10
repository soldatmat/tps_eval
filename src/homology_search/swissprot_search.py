from __future__ import annotations

"""Tool A — broad SEQUENCE homology search of designs vs Swiss-Prot.

For each design sequence, DIAMOND blastp against a DIAMOND DB built from the full
Swiss-Prot (`uniprot_sprot.fasta`), then report the top hit across ALL proteins and
whether it (and the top-N) are terpene synthases (TPS) — catching function-drift
(closest relative is a non-TPS enzyme) and confirming specificity.

DIAMOND is chosen over MMseqs2 here because: (a) it is a single self-contained
binary trivially installed via `conda install -c bioconda -c conda-forge diamond`,
(b) `diamond makedb` builds the Swiss-Prot DB in seconds from the same
`uniprot_sprot.fasta` we already need, and (c) blastp with `--very-sensitive` is
more than fast enough for the (few-thousand) design set on CPU. MMseqs2 would also
work; DIAMOND keeps the toolchain minimal and the accession parsing simple.

Output CSV keyed by ID (one row per design, NaN/empty if no hits):
    swissprot_top_hit            accession of the single best hit (max bitscore)
    swissprot_top_pident         % identity of that best hit
    swissprot_top_bitscore       bitscore of that best hit
    swissprot_top_is_tps         bool: is the best hit a TPS?
    swissprot_best_nontps_pident best pident among non-TPS hits (function-drift signal)
    swissprot_n_tps_in_topN      how many of the top-N hits are TPS
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers  # noqa: E402
from homology_search.tps_accessions import is_tps, load_tps_accessions  # noqa: E402

# Swiss-Prot FASTA headers are `sp|<ACCESSION>|<ENTRY_NAME> ...`. DIAMOND's default
# `sseqid` is the full `sp|ACC|NAME` token; we extract the middle field as the
# UniProt accession. `diamond makedb` keeps the header verbatim, so this holds.
DIAMOND_OUTFMT = ["qseqid", "sseqid", "pident", "bitscore", "evalue"]

COLUMNS = [
    "ID",
    "swissprot_top_hit",
    "swissprot_top_pident",
    "swissprot_top_bitscore",
    "swissprot_top_is_tps",
    "swissprot_best_nontps_pident",
    "swissprot_n_tps_in_topN",
]


def _accession_from_sseqid(sseqid: str) -> str:
    """Extract the UniProt accession from a DIAMOND subject id.

    Swiss-Prot: `sp|P12345|NAME_SPECIES` -> `P12345`. Falls back to the raw token
    if it is not pipe-delimited (e.g. a plain-accession DB)."""
    parts = str(sseqid).split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        return parts[1]
    return str(sseqid)


def _default_save_path(fasta_path: str) -> str:
    base, _ext = os.path.splitext(fasta_path)
    return f"{base}_swissprot_search.csv"


def _run_diamond(fasta_path: str, diamond_db: str, out_tsv: str, *,
                 top_n: int, threads: int, sensitivity: str) -> None:
    cmd = [
        "diamond", "blastp",
        "--db", diamond_db,
        "--query", fasta_path,
        "--out", out_tsv,
        "--outfmt", "6", *DIAMOND_OUTFMT,
        # Keep the top `top_n` hits per query (by score). --max-target-seqs caps it.
        "--max-target-seqs", str(top_n),
        f"--{sensitivity}",
        "--threads", str(threads),
        "--quiet",
    ]
    print("Running:", " ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end="")
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)


def _summarize_query(group: pd.DataFrame, tps_set: frozenset, top_n: int) -> dict:
    """Reduce one query's hit table (already top-N, sorted by bitscore desc) to a row."""
    group = group.sort_values("bitscore", ascending=False).head(top_n)
    accs = [_accession_from_sseqid(s) for s in group["sseqid"]]
    is_tps_flags = [is_tps(a, tps_set) for a in accs]

    top_acc = accs[0]
    top_is_tps = bool(is_tps_flags[0])

    nontps_pidents = [p for p, t in zip(group["pident"], is_tps_flags) if not t]
    best_nontps_pident = float(max(nontps_pidents)) if nontps_pidents else float("nan")

    return {
        "swissprot_top_hit": top_acc,
        "swissprot_top_pident": float(group["pident"].iloc[0]),
        "swissprot_top_bitscore": float(group["bitscore"].iloc[0]),
        "swissprot_top_is_tps": top_is_tps,
        "swissprot_best_nontps_pident": best_nontps_pident,
        "swissprot_n_tps_in_topN": int(sum(is_tps_flags)),
    }


def swissprot_search(
    fasta_path: str,
    diamond_db: str,
    tps_accessions_path: str,
    *,
    save_path: Optional[str] = None,
    top_n: int = 25,
    threads: int = 4,
    sensitivity: str = "very-sensitive",
) -> pd.DataFrame:
    """DIAMOND blastp every design vs Swiss-Prot; classify top hits TPS/non-TPS.

    Every design in `fasta_path` gets a row (NaN/empty if it had no hits)."""
    tps_set = load_tps_accessions(tps_accessions_path)
    identifiers, _seqs = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    # Match the other metrics' join key: ID is the first whitespace token.
    identifiers = [str(i).split(" ", 1)[0] for i in identifiers]

    tmp = tempfile.mkdtemp(prefix="swissprot_search_")
    out_tsv = os.path.join(tmp, "hits.tsv")
    try:
        _run_diamond(fasta_path, diamond_db, out_tsv,
                     top_n=top_n, threads=threads, sensitivity=sensitivity)
        if os.path.getsize(out_tsv) > 0:
            hits = pd.read_csv(out_tsv, sep="\t", header=None, names=DIAMOND_OUTFMT)
        else:
            hits = pd.DataFrame(columns=DIAMOND_OUTFMT)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # DIAMOND's qseqid is the design's FASTA id (first token) — align with `identifiers`.
    hits["qseqid"] = hits["qseqid"].astype(str).map(lambda x: x.split(" ", 1)[0])

    summaries = {}
    for qid, group in hits.groupby("qseqid"):
        summaries[qid] = _summarize_query(group, tps_set, top_n)

    rows: List[dict] = []
    empty = {c: (float("nan") if c not in ("swissprot_top_hit",) else "")
             for c in COLUMNS if c != "ID"}
    for ident in identifiers:
        row = {"ID": ident}
        row.update(summaries.get(ident, dict(empty)))
        # n_tps_in_topN is an int count; default 0 not NaN when no hits.
        if ident not in summaries:
            row["swissprot_n_tps_in_topN"] = 0
            row["swissprot_top_is_tps"] = pd.NA
        rows.append(row)

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    n_with_hits = int(df["swissprot_top_hit"].astype(str).str.len().gt(0).sum())
    print(f"Designs: {len(df)}  with >=1 hit: {n_with_hits}  "
          f"top-hit-is-TPS: {int((df['swissprot_top_is_tps'] == True).sum())}")

    if save_path is None:
        save_path = _default_save_path(fasta_path)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
