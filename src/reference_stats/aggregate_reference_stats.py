"""Aggregate per-metric reference CSVs into a single reference-stats JSON.

This is the back half of the *reference-statistics pipeline* (see
``scripts/compute_reference_stats.sh``). The applicable metric tools are run on
the MARTS-DB known-TPS reference set (their normal ``run_<tool>.sh`` runners,
which emit per-design CSVs keyed by ``ID``); this script reads those CSVs and
condenses each metric *column* into summary statistics describing the natural
TPS distribution. The result is a single committable JSON the rest of the
project loads to draw "natural TPS" bands around generated-design results.

Output shape (JSON)::

    {
      "<metric>": {                       # e.g. "motif_pair_distance"
        "n_rows": 1195,                   # rows in that metric's CSV
        "source_csv": "....csv",          # basename of the input CSV
        "columns": {
          "<column>": {                   # numeric column
            "kind": "numeric",
            "count": 1180, "n_missing": 15,
            "mean": ..., "std": ...,
            "min": ..., "p1": ..., "p5": ..., "p25": ...,
            "median": ..., "p75": ..., "p95": ..., "p99": ..., "max": ...
          },
          "<column>": {                   # categorical / boolean column
            "kind": "categorical",
            "count": 1195, "n_missing": 0, "n_unique": 12,
            "frequencies": {"alpha": 540, "alpha-beta": 320, "": 41, ...}
          }
        }
      },
      ...
    }

Design notes:
* **Column-type driven, not metric-hardcoded.** Any CSV keyed by ``ID`` is
  accepted; each non-ID column is classified at runtime as numeric (full stats)
  or categorical (frequency table). New metric tools (e.g. a freshly added
  ``radius_of_gyration``) are picked up automatically as long as their CSV is in
  the input dir — no code change here.
* **Helper / identifier columns are dropped.** ``ID`` is always dropped. String
  columns that are clearly per-design annotations rather than a distribution
  (matched-motif substrings, the chosen template name, etc.) are dropped too;
  see ``_DROP_COLUMNS``. Everything else that is non-numeric becomes a
  categorical frequency table (covers ``domain_architecture`` and any boolean
  motif-presence flags).
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# The percentile grid for numeric columns (besides min/max/mean/std/median).
_PERCENTILES = {
    "p1": 1,
    "p5": 5,
    "p25": 25,
    "p75": 75,
    "p95": 95,
    "p99": 99,
}

# Columns that are per-design annotations, NOT a distribution worth banding.
# Dropped from every CSV before classification. Matching is case-sensitive on the
# exact column name. ``ID`` is handled separately (always the key).
_DROP_COLUMNS = {
    "ID",
    # motif_pair_distance helper columns (matched substrings + raw start indices)
    "ddxxd_motif",
    "ddxxd_start",
    "nse_dte_motif",
    "nse_dte_start",
    # active_site_geometry: which catalytic template matched (a label, per design)
    "best_template",
    # prepare_csv / generic: raw sequence column if a tool ever carries it through
    "sequence",
}

# Map metric-name -> CSV filename suffix. The MARTS-DB reference CSVs are named
# ``<input>_<suffix>.csv`` by the tools (fasta-keyed for the sequence branch,
# structs-dir-keyed for the structure branch). We discover by suffix so the
# input dir can mix sequence- and structure-metric outputs.
METRIC_SUFFIXES: Dict[str, str] = {
    # sequence branch (intrinsic properties)
    "motif_pair_distance": "motif_pair_distance",
    "esm_pseudo_perplexity": "esm_pseudo_perplexity",
    # structure branch (intrinsic properties; computed only if structures exist)
    "plddt": "plddt",
    "motif_structural_distance": "motif_structural_distance",
    "active_site_geometry": "active_site_geometry",
    "aggregation": "aggregation",
    "domain_composition": "domain_composition",
    "proteinmpnn_score": "proteinmpnn_score",
    "radius_of_gyration": "radius_of_gyration",
}


def _numeric_stats(series: pd.Series) -> Dict[str, Optional[float]]:
    """Full summary stats for one numeric column (NaNs excluded from stats)."""
    values = pd.to_numeric(series, errors="coerce")
    n_total = int(len(values))
    valid = values.dropna()
    count = int(len(valid))
    stats: Dict[str, Optional[float]] = {
        "kind": "numeric",
        "count": count,
        "n_missing": n_total - count,
    }
    if count == 0:
        for key in ("mean", "std", "min", "median", "max", *_PERCENTILES):
            stats[key] = None
        return stats
    arr = valid.to_numpy(dtype=float)
    stats["mean"] = float(np.mean(arr))
    # ddof=1 sample std; undefined for a single point -> 0.0 by numpy convention
    stats["std"] = float(np.std(arr, ddof=1)) if count > 1 else 0.0
    stats["min"] = float(np.min(arr))
    stats["median"] = float(np.median(arr))
    stats["max"] = float(np.max(arr))
    pct_keys = list(_PERCENTILES)
    pct_vals = np.percentile(arr, [_PERCENTILES[k] for k in pct_keys])
    for key, val in zip(pct_keys, pct_vals):
        stats[key] = float(val)
    return stats


def _categorical_stats(series: pd.Series) -> Dict[str, object]:
    """Frequency table for one categorical / boolean column."""
    n_total = int(len(series))
    non_null = series.dropna()
    # Booleans / NaN-as-category: count empty-string and bool values verbatim.
    counts = non_null.astype(str).value_counts()
    freqs = {str(k): int(v) for k, v in counts.items()}
    return {
        "kind": "categorical",
        "count": int(len(non_null)),
        "n_missing": n_total - int(len(non_null)),
        "n_unique": int(len(freqs)),
        "frequencies": freqs,
    }


def _classify_and_stat(series: pd.Series) -> Dict[str, object]:
    """Decide numeric vs categorical for one column and compute its stats."""
    # A column is numeric if pandas already typed it numeric, OR if coercion
    # leaves at least one real number and no genuine string labels. Object
    # columns holding category labels (domain_architecture, motif flags) coerce
    # entirely to NaN -> treated as categorical.
    if pd.api.types.is_bool_dtype(series):
        return _categorical_stats(series)
    if pd.api.types.is_numeric_dtype(series):
        return _numeric_stats(series)
    coerced = pd.to_numeric(series, errors="coerce")
    # If everything that is non-null coerces to a number, it's a numeric column
    # stored as object (rare). Otherwise treat as categorical.
    non_null = series.dropna()
    if len(non_null) > 0 and coerced.dropna().shape[0] == non_null.shape[0]:
        return _numeric_stats(coerced)
    return _categorical_stats(series)


def aggregate_csv(csv_path: str) -> Dict[str, object]:
    """Aggregate one metric CSV into ``{n_rows, source_csv, columns: {...}}``."""
    df = pd.read_csv(csv_path)
    columns: Dict[str, object] = {}
    for col in df.columns:
        if col in _DROP_COLUMNS:
            continue
        columns[col] = _classify_and_stat(df[col])
    return {
        "n_rows": int(len(df)),
        "source_csv": os.path.basename(csv_path),
        "columns": columns,
    }


def discover_csvs(input_dir: str) -> Dict[str, str]:
    """Map metric-name -> CSV path for every known metric present in input_dir.

    Matches ``*_<suffix>.csv`` (the tools' ``<input>_<suffix>.csv`` naming). If
    several files match a suffix, the lexicographically last is used and a note
    is printed. Unknown ``*.csv`` files in the dir are ignored.
    """
    found: Dict[str, str] = {}
    for metric, suffix in METRIC_SUFFIXES.items():
        matches = sorted(glob.glob(os.path.join(input_dir, f"*_{suffix}.csv")))
        # exact "<suffix>.csv" (no leading underscore) also accepted
        matches += sorted(glob.glob(os.path.join(input_dir, f"{suffix}.csv")))
        matches = sorted(set(matches))
        if not matches:
            continue
        if len(matches) > 1:
            print(f"[warn] {metric}: multiple CSVs match, using last: {matches}")
        found[metric] = matches[-1]
    return found


def build_reference_stats(
    input_dir: str,
    *,
    explicit: Optional[Dict[str, str]] = None,
) -> Dict[str, object]:
    """Build the full reference-stats dict from a dir of per-metric CSVs.

    ``explicit`` optionally maps metric-name -> csv path, overriding/augmenting
    suffix-based discovery (used when a CSV lives elsewhere or is renamed).
    """
    csvs = discover_csvs(input_dir)
    if explicit:
        csvs.update(explicit)
    if not csvs:
        raise SystemExit(
            f"No recognised metric CSVs found in {input_dir!r}. "
            f"Expected files like *_motif_pair_distance.csv. "
            f"Known suffixes: {sorted(METRIC_SUFFIXES.values())}"
        )
    out: Dict[str, object] = {}
    for metric in sorted(csvs):
        path = csvs[metric]
        print(f"[ok] {metric}: {path}")
        out[metric] = aggregate_csv(path)
    return out


def _json_safe(obj):
    """Replace NaN/Inf with None so the JSON is strict-valid."""
    if isinstance(obj, float):
        return None if (math.isnan(obj) or math.isinf(obj)) else obj
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_safe(v) for v in obj]
    return obj


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate per-metric reference CSVs (computed on the "
        "MARTS-DB known-TPS set) into a single reference-stats JSON keyed by "
        "metric -> column -> stats. Numeric columns get full distribution stats "
        "(count, mean, std, min, p1/p5/p25, median, p75/p95/p99, max); "
        "categorical columns (e.g. domain_architecture) get a frequency table."
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing the per-metric reference CSVs "
        "(e.g. *_motif_pair_distance.csv).",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "marts_db_metric_stats.json",
        ),
        help="Output JSON path (default: src/reference_stats/marts_db_metric_stats.json).",
    )
    parser.add_argument(
        "--reference_name",
        default="marts_db",
        help="Reference set label embedded in the JSON metadata (default: marts_db).",
    )
    args = parser.parse_args()

    stats = build_reference_stats(args.input_dir)
    document = {
        "reference_set": args.reference_name,
        "input_dir": os.path.abspath(args.input_dir),
        "metrics": _json_safe(stats),
    }
    with open(args.output, "w") as fh:
        json.dump(document, fh, indent=2, sort_keys=False)
        fh.write("\n")
    n_metrics = len(stats)
    print(f"\nWrote {n_metrics} metric(s) to {args.output}")


if __name__ == "__main__":
    main()
