from __future__ import annotations

"""Tool B — broad STRUCTURE homology search of designs vs AlphaFold-Swiss-Prot.

For each design structure, foldseek easy-search against the ANNOTATED AlphaFold DB
Swiss-Prot set (`Alphafold/Swiss-Prot`), then report the top hit across ALL proteins
(by TM-score) and whether it (and the top-N) are terpene synthases — the structural
analog of the sequence search. Catches function-drift on the fold level and confirms
the design's structure is specifically TPS-like rather than a related-but-different
fold (e.g. a prenyltransferase).

The foldseek AFDB-Swiss-Prot target IDs encode the UniProt accession as
`AF-<ACC>-F1-model_v4`; we extract `<ACC>` and classify via the shared TPS set.

Output CSV keyed by ID (filename stem / af_output job name; reuses plddt.py's
af3-vs-flat detection and dir-keyed naming). One row per design, NaN/empty if no hits:
    foldseek_sprot_top_hit            accession of the best hit (max alntmscore)
    foldseek_sprot_top_tmscore        alntmscore of that best hit
    foldseek_sprot_top_is_tps         bool: is the best hit a TPS?
    foldseek_sprot_best_nontps_tmscore best alntmscore among non-TPS hits (drift signal)
    foldseek_sprot_n_tps_in_topN      how many of the top-N hits are TPS
"""

import os
import re
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

# Reuse the af3-vs-flat input detection + dir-keyed naming from the pLDDT tool.
from structure_metrics.plddt import _collect_structures  # noqa: E402
from homology_search.tps_accessions import is_tps, load_tps_accessions  # noqa: E402

# AFDB target ids look like `AF-<ACC>-F<frag>-model_v<n>` (sometimes with a `.pdb`
# / `.cif` suffix). Capture the accession.
_AFDB_RE = re.compile(r"^AF-([0-9A-Za-z]+)-F\d+-model")

FOLDSEEK_OUTFMT = "query,target,fident,alnlen,evalue,bits,alntmscore"
FOLDSEEK_COLS = ["query", "target", "fident", "alnlen", "evalue", "bits", "alntmscore"]

COLUMNS = [
    "ID",
    "foldseek_sprot_top_hit",
    "foldseek_sprot_top_tmscore",
    "foldseek_sprot_top_is_tps",
    "foldseek_sprot_best_nontps_tmscore",
    "foldseek_sprot_n_tps_in_topN",
]


def _accession_from_target(target: str) -> str:
    """Extract the UniProt accession from an AFDB foldseek target id.

    `AF-P12345-F1-model_v4` -> `P12345`. Falls back to the raw token (minus any
    structure-file extension) if it doesn't match the AFDB pattern."""
    t = str(target)
    m = _AFDB_RE.match(t)
    if m:
        return m.group(1)
    # Strip a trailing extension if present, else return verbatim.
    return os.path.splitext(t)[0]


def _query_to_id(query: str) -> str:
    """foldseek's `query` is the input structure filename; map it to the design ID
    (the filename stem, matching _collect_structures' keys and plddt.py)."""
    stem = os.path.basename(str(query))
    # Strip known structure extensions (handles `.pdb`, `.cif`, `.mmcif`, plus the
    # AF3 `_model.cif` whose stem we want to be the job name).
    for ext in (".pdb", ".cif", ".mmcif"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    if stem.endswith("_model"):
        stem = stem[: -len("_model")]
    return stem


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_foldseek_swissprot_search.csv")


def _run_foldseek(query_dir: str, afdb_db: str, out_tsv: str, tmp_dir: str,
                  *, max_seqs: int) -> None:
    cmd = (
        f"foldseek easy-search {query_dir} {afdb_db} {out_tsv} {tmp_dir} "
        f"--max-seqs {max_seqs} --format-output {FOLDSEEK_OUTFMT} -v 3"
    ).split()
    print("Running:", " ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end="")
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)


def _summarize_query(group: pd.DataFrame, tps_set: frozenset, top_n: int) -> dict:
    """Reduce one query's hit table to a row, ranked by alntmscore desc."""
    group = group.sort_values("alntmscore", ascending=False).head(top_n)
    accs = [_accession_from_target(t) for t in group["target"]]
    is_tps_flags = [is_tps(a, tps_set) for a in accs]

    top_acc = accs[0]
    nontps_tm = [tm for tm, t in zip(group["alntmscore"], is_tps_flags) if not t]
    best_nontps_tm = float(max(nontps_tm)) if nontps_tm else float("nan")

    return {
        "foldseek_sprot_top_hit": top_acc,
        "foldseek_sprot_top_tmscore": float(group["alntmscore"].iloc[0]),
        "foldseek_sprot_top_is_tps": bool(is_tps_flags[0]),
        "foldseek_sprot_best_nontps_tmscore": best_nontps_tm,
        "foldseek_sprot_n_tps_in_topN": int(sum(is_tps_flags)),
    }


def foldseek_swissprot_search(
    structs_dir: str,
    afdb_db: str,
    tps_accessions_path: str,
    *,
    save_path: Optional[str] = None,
    top_n: int = 25,
    max_seqs: int = 300,
) -> pd.DataFrame:
    """foldseek every design structure vs AFDB-Swiss-Prot; classify top hits TPS/non-TPS.

    Every design (af3 job name or flat-dir stem) gets a row, NaN/empty if no hits."""
    tps_set = load_tps_accessions(tps_accessions_path)
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    tmp = tempfile.mkdtemp(prefix="foldseek_sprot_")
    out_tsv = os.path.join(tmp, "hits.tsv")
    fs_tmp = os.path.join(tmp, "fs_tmp")
    # In af3 mode the structures live in per-job subfolders; foldseek can't take a
    # nested tree, so stage the chosen model files into one flat query dir, named by ID.
    query_dir = os.path.join(tmp, "query")
    os.makedirs(query_dir, exist_ok=True)
    try:
        for ident, path in structures.items():
            ext = os.path.splitext(path)[1] or ".pdb"
            shutil.copy(path, os.path.join(query_dir, f"{ident}{ext}"))

        _run_foldseek(query_dir, afdb_db, out_tsv, fs_tmp, max_seqs=max_seqs)
        if os.path.isfile(out_tsv) and os.path.getsize(out_tsv) > 0:
            hits = pd.read_csv(out_tsv, sep="\t", header=None, names=FOLDSEEK_COLS)
        else:
            hits = pd.DataFrame(columns=FOLDSEEK_COLS)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    hits["_id"] = hits["query"].map(_query_to_id)
    summaries = {}
    for qid, group in hits.groupby("_id"):
        summaries[qid] = _summarize_query(group, tps_set, top_n)

    rows: List[dict] = []
    for ident in structures.keys():
        row = {"ID": ident}
        if ident in summaries:
            row.update(summaries[ident])
        else:
            row.update({
                "foldseek_sprot_top_hit": "",
                "foldseek_sprot_top_tmscore": float("nan"),
                "foldseek_sprot_top_is_tps": pd.NA,
                "foldseek_sprot_best_nontps_tmscore": float("nan"),
                "foldseek_sprot_n_tps_in_topN": 0,
            })
        rows.append(row)

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    n_with_hits = int(df["foldseek_sprot_top_hit"].astype(str).str.len().gt(0).sum())
    print(f"Designs: {len(df)}  with >=1 hit: {n_with_hits}  "
          f"top-hit-is-TPS: {int((df['foldseek_sprot_top_is_tps'] == True).sum())}")

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
