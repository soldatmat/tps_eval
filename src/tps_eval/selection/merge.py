"""Merge the pipeline's per-tool CSVs into one wide table keyed by ID.

The selection layer (gate / band_filter / score / diversity_dedup / select) operates
on ONE wide table of every metric a design has, keyed by ID. This module produces that
table from the scattered ``<base>_<tool>.csv`` (sequence branch) and
``<structs_dir>_<tool>.csv`` (structure branch) outputs the eval pipeline writes.

It intentionally mirrors the merge conventions in ``src/dashboard/build_dashboard.py``
(``load_design_batch``) — same ID-column candidates, same ``_self`` column-suffixing to
keep gen-vs-self distinct from gen-vs-train, same missing-value tokens. The two differ
only in the target shape: the dashboard wants dict-of-column-lists for its SVG overlay;
selection wants a pandas DataFrame it can filter/rank. Keep the two conventions in sync
(if you add a metric whose CSV shares a column name with another tool, add its ``_self``
handling here too).

Unlike the dashboard, this KEEPS the ``sequence`` column (selection needs it to emit the
survivor FASTA and to feed diversity_dedup / order-preparation) and never subsamples.
"""
from __future__ import annotations

import glob
import os
from typing import List, Optional

import pandas as pd

# ID column candidates, in priority order (mirrors build_dashboard). The row-id column
# name varies by tool: ID / lowercase id / fa_id (SoluProt) / reference_id.
ID_CANDIDATES = ("ID", "id", "Id", "fa_id", "reference_id")

# Values that mean "missing" across the tools' CSVs (mirrors build_dashboard._MISSING).
MISSING_TOKENS = ("", "nan", "NA", "NaN", "None", "null")

# Pure identifier / duplicate-key columns dropped from every non-first file so they don't
# collide. NOTE: unlike the dashboard we KEEP `sequence` (needed to emit FASTAs) — it is
# de-duplicated by the first-wins column rule below, not dropped.
_ID_LIKE = set(ID_CANDIDATES)

# A CSV with more data columns than this is a raw feature matrix (e.g. the 1280-dim ESM
# embedding), not a per-design metric table — skip it (mirrors build_dashboard).
_MAX_COLUMNS_PER_FILE = 256


def resolve_csv_paths(entries: List[str]) -> List[str]:
    """Expand a list of files / directories / globs into a deduped, ordered CSV list.

    A directory contributes its ``*.csv`` children; a glob its matches; a file itself.
    Order is preserved (first occurrence wins on merge), duplicates removed.
    """
    paths: List[str] = []
    seen = set()
    for entry in entries:
        matches: List[str] = []
        if os.path.isdir(entry):
            matches = sorted(glob.glob(os.path.join(entry, "*.csv")))
        elif any(c in entry for c in "*?[]"):
            matches = sorted(glob.glob(entry))
        elif os.path.isfile(entry):
            matches = [entry]
        for p in matches:
            ap = os.path.abspath(p)
            if ap not in seen and p.endswith(".csv"):
                seen.add(ap)
                paths.append(p)
    return paths


def _id_column(columns) -> Optional[str]:
    return next((c for c in ID_CANDIDATES if c in columns), None)


def merge_metrics(entries: List[str]) -> pd.DataFrame:
    """Merge the per-tool CSVs named by ``entries`` into one wide DataFrame indexed by ID.

    - Each CSV's ID column (first of ID_CANDIDATES present) becomes the join key.
    - Self-comparison tools (``*_self.csv``) get their data columns suffixed ``_self`` so
      they don't clobber their gen-vs-train sibling (mirrors the dashboard convention).
    - Collisions are resolved per-(ID, column) FIRST-WINS over the UNION of IDs (mirrors
      the dashboard's ``load_design_batch``): a column name repeated across files — the
      same tool run on several input FASTAs with DISJOINT IDs (e.g. one per class/setting),
      or two tools sharing ``sequence`` — is unioned by ID, each cell taken from the first
      file that has a value for it. (A naive column-level first-wins would blank out the
      rows only present in the later files.)
    - Missing-value tokens are normalised to NaN; numeric-looking columns are coerced to
      float, the rest left as strings.

    Raises ValueError if no usable CSVs are found.
    """
    csv_paths = resolve_csv_paths(entries)
    if not csv_paths:
        raise ValueError(f"no usable CSVs found in: {entries}")

    merged: Optional[pd.DataFrame] = None
    for path in csv_paths:
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False)
        except (pd.errors.EmptyDataError, pd.errors.ParserError):
            continue
        id_col = _id_column(df.columns)
        if id_col is None or df.empty:
            continue
        data_cols = [c for c in df.columns if c != id_col and c not in _ID_LIKE]
        if len(data_cols) > _MAX_COLUMNS_PER_FILE:
            print(f"  [merge] skipping {os.path.basename(path)}: {len(data_cols)} data "
                  f"columns (> {_MAX_COLUMNS_PER_FILE}) — looks like a raw feature matrix.")
            continue
        # Self-file columns get a `_self` suffix (unless already so named).
        tool_is_self = os.path.basename(path).rsplit(".", 1)[0].endswith("_self")
        rename = {c: (f"{c}_self" if tool_is_self and not c.endswith("_self") else c)
                  for c in data_cols}
        sub = df[[id_col] + data_cols].rename(columns={**{id_col: "ID"}, **rename})
        sub = sub.drop_duplicates(subset="ID").set_index("ID")
        # Normalise missing tokens to NaN BEFORE combining, else an earlier file's blank
        # (kept-default-na off -> "") would count as present and win over a later real value.
        sub = sub.replace(list(MISSING_TOKENS), pd.NA)
        # combine_first: keep existing (earlier-file) values per (ID, col), fill gaps from
        # this file, and UNION the index + columns. -> per-cell first-wins across all files.
        merged = sub if merged is None else merged.combine_first(sub)

    if merged is None:
        raise ValueError(f"no CSV under {entries} had a recognised ID column {ID_CANDIDATES}.")

    # Coerce numeric-looking columns to float (the rest stay strings).
    for c in merged.columns:
        coerced = pd.to_numeric(merged[c], errors="coerce")
        # Treat a column as numeric only if every non-missing value parsed (else keep str).
        nonnull = merged[c].notna()
        if nonnull.any() and coerced[nonnull].notna().all():
            merged[c] = coerced
    return merged.reset_index()


def write_merged(df: pd.DataFrame, output_path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"[merge] wrote {len(df)} rows x {len(df.columns)} cols -> {output_path}")
