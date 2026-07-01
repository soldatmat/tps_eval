"""diversity_dedup — reduce to a diverse, quality-prioritised subset.

Clusters the designs by sequence identity (MMseqs2 ``easy-cluster``) at a chosen threshold,
keeps the BEST-quality representative of each cluster, then takes the top-N by quality. This
is diversity-CONSTRAINED, quality-PRIORITISED selection (not maximal spread): every survivor
is the best member of its own identity cluster, and no two survivors exceed the clustering
identity. The threshold may be set PER GROUP (e.g. per first-cyclization class), since some
classes are more diversity-limited than others.

MMseqs2 is the same tool the local_sequence_search metric already uses; ``easy-cluster`` is
invoked per group. Higher ``quality_col`` is better.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from typing import Dict, Optional, Tuple

import pandas as pd

from io_fasta import write_fasta


def _mmseqs_clusters(id_to_seq: Dict[str, str], min_seq_id: float,
                     coverage: float = 0.8, threads: int = 4) -> Dict[str, str]:
    """Return {member_id: cluster_representative_id} from MMseqs2 easy-cluster."""
    with tempfile.TemporaryDirectory() as work:
        fasta = os.path.join(work, "seqs.fasta")
        write_fasta(id_to_seq, fasta)
        out_prefix = os.path.join(work, "clu")
        tmp = os.path.join(work, "tmp")
        cmd = ["mmseqs", "easy-cluster", fasta, out_prefix, tmp,
               "--min-seq-id", str(min_seq_id), "-c", str(coverage),
               "--threads", str(threads)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"mmseqs easy-cluster failed:\n{proc.stdout}\n{proc.stderr}")
        tsv = out_prefix + "_cluster.tsv"  # cols: representative<TAB>member
        member_to_rep: Dict[str, str] = {}
        with open(tsv) as fh:
            for line in fh:
                rep, member = line.rstrip("\n").split("\t")[:2]
                member_to_rep[member.split()[0]] = rep.split()[0]
        return member_to_rep


def _dedup_group(g: pd.DataFrame, id_col: str, seq_col: str, quality_col: str,
                 min_seq_id: float, n_out: Optional[int], coverage: float) -> pd.DataFrame:
    id_to_seq = {str(r[id_col]): str(r[seq_col]) for _, r in g.iterrows()
                 if pd.notna(r[seq_col])}
    if len(id_to_seq) <= 1:
        reps = g
    else:
        member_to_rep = _mmseqs_clusters(id_to_seq, min_seq_id, coverage)
        g = g.copy()
        g["_cluster"] = g[id_col].astype(str).map(member_to_rep).fillna(g[id_col].astype(str))
        # Best-quality representative per cluster.
        reps = (g.sort_values(quality_col, ascending=False)
                .drop_duplicates(subset="_cluster", keep="first")
                .drop(columns=["_cluster"]))
    reps = reps.sort_values(quality_col, ascending=False)
    if n_out is not None:
        reps = reps.head(n_out)
    return reps


def apply_diversity_dedup(df: pd.DataFrame, quality_col: str,
                          id_threshold: Optional[float] = None,
                          id_threshold_per_group: Optional[Dict[str, float]] = None,
                          group_col: Optional[str] = None,
                          n_out_per_group: Optional[int] = None,
                          seq_col: str = "sequence", id_col: str = "ID",
                          coverage: float = 0.8) -> Tuple[pd.DataFrame, Dict]:
    """Cluster-dedup per group at the group's identity threshold, keep best-rep, top-N."""
    if seq_col not in df.columns:
        raise ValueError(f"diversity_dedup needs a '{seq_col}' column (got {list(df.columns)}).")
    if quality_col not in df.columns:
        raise ValueError(f"diversity_dedup quality_col '{quality_col}' not present.")

    def thr_for(group_key) -> float:
        if id_threshold_per_group and group_key is not None:
            t = id_threshold_per_group.get(str(group_key))
            if t is not None:
                return t
        if id_threshold is None:
            raise ValueError("diversity_dedup needs id_threshold or id_threshold_per_group "
                             f"(no threshold for group '{group_key}').")
        return id_threshold

    per_group = []
    if group_col and group_col in df.columns:
        parts = []
        for key, g in df.groupby(group_col, sort=False):
            kept = _dedup_group(g, id_col, seq_col, quality_col, thr_for(key),
                                n_out_per_group, coverage)
            per_group.append({"group": str(key), "n_in": len(g), "n_out": len(kept),
                              "min_seq_id": thr_for(key)})
            parts.append(kept)
        out = pd.concat(parts) if parts else df.iloc[0:0]
    else:
        out = _dedup_group(df, id_col, seq_col, quality_col, thr_for(None),
                           n_out_per_group, coverage)
        per_group.append({"group": "all", "n_in": len(df), "n_out": len(out),
                          "min_seq_id": thr_for(None)})
    report = {"op": "diversity_dedup", "n_in": len(df), "n_out": len(out),
              "quality_col": quality_col, "groups": per_group}
    return out, report
