"""select — the composite selection driver.

Runs an ordered pipeline of selection ops (gate / band_filter / score / diversity_dedup)
described by a JSON spec over a merged metric table, then takes the top-N per group. Emits:
  * ``<prefix>_survivors.csv``   — the surviving rows (all merged columns + any score),
  * ``<prefix>_survivors.fasta`` — their sequences (for the next funnel tier / order prep),
  * ``<prefix>_manifest.md``     — a provenance record: every op, its parameters, and the
    in/out counts per group (auto-generates the SELECTION_PROCEDURE-style record we used
    to write by hand).

Spec schema::

    {
      "group_by": "class",              # grouping column (optional)
      "group_from_id": "_(c\\d+)_",     # optional regex to synthesise group_by from the ID
      "n_out_per_group": 16,            # final per-group cap (by score if present)
      "ops": [ {op-spec}, ... ]         # run in order; see each primitive's docstring
    }

Op specs:
  gate:            {"op":"gate", "conditions":[...]}
  band_filter:     {"op":"band_filter", "metrics":{...}, "bands_file":"..."}
  score:           {"op":"score", "terms":[...], "zscore_within":"class"}   # default group_by
  diversity_dedup: {"op":"diversity_dedup", "quality_col":"score",
                    "id_threshold":0.7 | "id_threshold_per_group":{...}}
"""
from __future__ import annotations

import json
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from band_filter import apply_band_filter
from diversity_dedup import apply_diversity_dedup
from gate import apply_gate
from io_fasta import read_fasta_map, write_fasta
from score import apply_score

_SCORE_COL = "score"


def _synthesise_group(df: pd.DataFrame, spec: dict) -> Tuple[pd.DataFrame, Optional[str]]:
    group_by = spec.get("group_by")
    regex = spec.get("group_from_id")
    if group_by and regex and group_by not in df.columns:
        pat = re.compile(regex)
        def _extract(rid):
            m = pat.search(str(rid))
            return m.group(1) if m else None
        df = df.copy()
        df[group_by] = df["ID"].map(_extract)
        missing = int(df[group_by].isna().sum())
        if missing:
            print(f"  [select] group_from_id '{regex}' matched no group for {missing} IDs.")
    return df, group_by


def run_selection(df: pd.DataFrame, spec: dict,
                  fasta_map: Optional[Dict[str, str]] = None
                  ) -> Tuple[pd.DataFrame, List[dict], Optional[str]]:
    """Apply the spec's ops in order, then cap to n_out_per_group. Returns
    (survivors_df, op_reports, group_col)."""
    df, group_by = _synthesise_group(df, spec)
    # Inject sequences up front (from the seed FASTA) so ops that need them mid-pipeline —
    # notably diversity_dedup — can use them, not just the final FASTA write.
    if fasta_map is not None and "sequence" not in df.columns:
        df = df.copy()
        df["sequence"] = df["ID"].map(fasta_map)
    reports: List[dict] = [{"op": "input", "n_in": len(df),
                            "group_counts": _group_counts(df, group_by)}]
    cur = df
    for op_spec in spec.get("ops", []):
        op = op_spec["op"]
        if op == "gate":
            cur, rep = apply_gate(cur, op_spec["conditions"])
        elif op == "band_filter":
            cur, rep = apply_band_filter(cur, op_spec.get("metrics", {}),
                                         bands_file=op_spec.get("bands_file"))
        elif op == "score":
            cur, rep = apply_score(cur, op_spec["terms"],
                                   zscore_within=op_spec.get("zscore_within", group_by),
                                   score_col=_SCORE_COL)
        elif op == "diversity_dedup":
            cur, rep = apply_diversity_dedup(
                cur, quality_col=op_spec.get("quality_col", _SCORE_COL),
                id_threshold=op_spec.get("id_threshold"),
                id_threshold_per_group=op_spec.get("id_threshold_per_group"),
                group_col=op_spec.get("group_col", group_by),
                n_out_per_group=op_spec.get("n_out_per_group", spec.get("n_out_per_group")),
                seq_col=op_spec.get("seq_col", "sequence"),
                coverage=op_spec.get("coverage", 0.8))
        else:
            raise ValueError(f"unknown selection op '{op}'")
        rep["group_counts"] = _group_counts(cur, group_by)
        reports.append(rep)

    # Final per-group cap by score (if a score column exists) — no-op when diversity_dedup
    # already produced <= n_out per group, so the two compose idempotently.
    n_out = spec.get("n_out_per_group")
    if n_out is not None:
        cur = _take_top_n(cur, group_by, n_out)
        reports.append({"op": "take_top_n", "n_out_per_group": n_out, "n_out": len(cur),
                        "group_counts": _group_counts(cur, group_by)})

    if fasta_map is not None and "sequence" not in cur.columns:
        cur = cur.copy()
        cur["sequence"] = cur["ID"].map(fasta_map)
    return cur.reset_index(drop=True), reports, group_by


def _take_top_n(df: pd.DataFrame, group_by: Optional[str], n_out: int) -> pd.DataFrame:
    sort_col = _SCORE_COL if _SCORE_COL in df.columns else None
    if group_by and group_by in df.columns:
        if sort_col:
            df = df.sort_values(sort_col, ascending=False)
        return df.groupby(group_by, sort=False, group_keys=False).head(n_out)
    if sort_col:
        df = df.sort_values(sort_col, ascending=False)
    return df.head(n_out)


def _group_counts(df: pd.DataFrame, group_by: Optional[str]) -> Dict[str, int]:
    if group_by and group_by in df.columns:
        return {str(k): int(v) for k, v in df[group_by].value_counts(dropna=False).items()}
    return {"all": len(df)}


def _fmt_groups(groups: Dict[str, int]) -> str:
    return ", ".join(f"{k}={v}" for k, v in sorted(groups.items()))


def write_manifest(reports: List[dict], spec: dict, output_path: str,
                   title: str = "Selection") -> None:
    lines = [f"# {title} — provenance manifest", "",
             "Auto-generated by `src/selection/select_designs.py`. Each row is one selection "
             "op with its surviving count (per group).", "",
             "## Spec", "```json", json.dumps(spec, indent=2), "```", "",
             "## Funnel", "", "| step | op | survivors | per-group |",
             "|---|---|---|---|"]
    for r in reports:
        op = r["op"]
        n = r.get("n_out", r.get("n_pass", r.get("n_in", "")))
        detail = ""
        if op == "gate":
            detail = "; ".join(f"{c['condition']} → {c['passed']}" for c in r.get("conditions", []))
        elif op == "band_filter":
            detail = "; ".join(f"{m['metric']} → {m['passed']}" for m in r.get("metrics", []))
        elif op == "score":
            detail = "score = " + " + ".join(
                f"{'−' if t['direction']=='lower' else ''}z({t['col']})×{t['weight']}"
                for t in r.get("terms", []))
        elif op == "diversity_dedup":
            detail = "; ".join(f"{g['group']} @id≤{g['min_seq_id']}: {g['n_in']}→{g['n_out']}"
                               for g in r.get("groups", []) if isinstance(g, dict))
        lines.append(f"| {op} | {op} | {n} | {_fmt_groups(r.get('group_counts', {}))} |")
        if detail:
            lines.append(f"| | ↳ | | {detail} |")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[select] wrote manifest -> {output_path}")


def select_and_write(df: pd.DataFrame, spec: dict, output_prefix: str,
                     fasta_map: Optional[Dict[str, str]] = None,
                     title: str = "Selection") -> pd.DataFrame:
    survivors, reports, group_by = run_selection(df, spec, fasta_map=fasta_map)
    csv_path = output_prefix + "_survivors.csv"
    survivors.to_csv(csv_path, index=False)
    print(f"[select] {reports[0]['n_in']} → {len(survivors)} survivors "
          f"({_fmt_groups(_group_counts(survivors, group_by))}) -> {csv_path}")
    if "sequence" in survivors.columns:
        id_to_seq = {str(r["ID"]): str(r["sequence"]) for _, r in survivors.iterrows()
                     if pd.notna(r["sequence"])}
        if id_to_seq:
            write_fasta(id_to_seq, output_prefix + "_survivors.fasta")
            print(f"[select] wrote {len(id_to_seq)} sequences -> "
                  f"{output_prefix}_survivors.fasta")
    else:
        print("[select] no 'sequence' column and no --fasta given -> FASTA not written.")
    write_manifest(reports, spec, output_prefix + "_manifest.md", title=title)
    return survivors
