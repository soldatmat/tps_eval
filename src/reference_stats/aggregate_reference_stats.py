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
            "median": ..., "p75": ..., "p95": ..., "p99": ..., "max": ...,
            "by_<labeling>": {            # only when --group_by given
              "<class>": { ... same stats, restricted to rows in that class ... },
              ...
            }
          },
          "<column>": {                   # categorical / boolean column
            "kind": "categorical",
            "count": 1195, "n_missing": 0, "n_unique": 12,
            "frequencies": {"alpha": 540, "alpha-beta": 320, "": 41, ...},
            "by_<labeling>": { "<class>": { ... per-class freq table ... }, ... }
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
* **Per-class stratification (``--group_by``).** Optionally pass one or more
  ``reference_id,label`` CSV label files. For every metric column, the same
  stats are *also* computed restricted to each class, emitted under a
  ``by_<labeling>`` block keyed by class label. The overall (ungrouped) stats
  are kept alongside. This is label-agnostic and metric-agnostic: the labeling
  name comes from the file basename (overridable), the join is a plain
  ``ID``-vs-``reference_id`` string match, and a labeling whose IDs come from
  another metric's CSV (e.g. ``domain_architecture`` from
  ``*_domain_composition.csv``) is supported via the same mechanism.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
from typing import Dict, List, Optional, Tuple

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
    # row-index / leftover id-alias columns (after _normalize_id_column picks the
    # canonical ID, any other identifier column is a per-design annotation, not a band)
    "runtime_id",
    "fa_id",
    "id",
    "reference_id",
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

# Canonical metric-name -> CSV filename suffix, used only for PRETTY NAMING of the
# known metrics. Discovery itself is generic (see discover_csvs): every ``*.csv``
# keyed by ID in the input dir is banded EXCEPT labeling files and the inherently
# comparative metrics below. The MARTS-DB reference CSVs are named
# ``<input>_<suffix>.csv`` by the tools (fasta-keyed for the sequence branch,
# structs-dir-keyed for the structure branch).
METRIC_SUFFIXES: Dict[str, str] = {
    # sequence branch (intrinsic properties)
    "motif_pair_distance": "motif_pair_distance",
    "esm_pseudo_perplexity": "esm_pseudo_perplexity",
    "motif_search": "motifs",
    "soluprot": "soluprot",
    "enzyme_explorer_sequence_only": "enzyme_explorer_sequence_only",
    # structure branch (intrinsic properties; computed only if structures exist)
    "plddt": "plddt",
    "motif_structural_distance": "motif_structural_distance",
    "active_site_geometry": "active_site_geometry",
    "aggregation": "aggregation",
    "domain_composition": "domain_composition",
    "proteinmpnn_score": "proteinmpnn_score",
    "radius_of_gyration": "radius_of_gyration",
    "pocket_descriptors": "pocket_descriptors",
    "aromatic_lining": "aromatic_lining",
    "diphosphate_sensor": "diphosphate_sensor",
    # PAE-derived fold-confidence (intrinsic; needs a saved PAE npz)
    "global_confidence": "global_confidence",
    "interdomain_pae": "interdomain_pae",
    # active-site ion placement (only carries signal on holo folds with modelled
    # ions -- AF3 holo / Boltz2 holo; apo/ESMFold rows are not-applicable but
    # harmless to band)
    "ion_site_check": "ion_site_check",
    # prenyl-PP substrate positioning (holo-only: requires a co-folded substrate,
    # e.g. AF3 co-fold or Boltz2 holo; apo folds yield all-NaN, not-applicable)
    "substrate_positioning": "substrate_positioning",
}

# Inherently COMPARATIVE metrics: each measures similarity to a SEPARATE reference
# set (a train set, the known-TPS set, or Swiss-Prot), or is a leave-one-out
# prediction/transfer. A "natural band" is ill-defined for them on the reference set
# itself (the reference set IS the natural set), so they are excluded from banding
# even if a CSV is present. Matched as a metric name OR a trailing ``_<suffix>``.
_COMPARATIVE_SUFFIXES = {
    "max_sequence_identity", "max_sequence_identity_self",
    "local_sequence_search", "local_sequence_search_self",
    "embedding_esm1b", "embedding_esm1b_min_embedding_distance",
    "embedding_esm1b_min_embedding_distance_self",
    "min_embedding_distance", "min_embedding_distance_self",
    "swissprot_search", "foldseek_swissprot_search",
    "structural_identity", "domain_structural_identity",
    "knn_label_transfer", "sdr_divergence", "substrate_class",
    "self_consistency",
}

# Input-name prefixes the tools prepend; stripped to name an UNKNOWN (newly-added)
# metric's CSV when it isn't in METRIC_SUFFIXES, so new tools band with no code change.
_KNOWN_PREFIXES = (
    "structs_esmfold_", "structs_af3_", "structs_boltz2_holo_", "structs_boltz2_",
    "structs_", "TPS_sequences_",
)


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


def _column_stats(
    series: pd.Series,
    groupings: Optional[Dict[str, pd.Series]] = None,
) -> Dict[str, object]:
    """Stats for one column (overall) plus optional per-class ``by_<labeling>``.

    ``groupings`` maps ``labeling_name -> label Series`` aligned to ``series``'s
    index (same row order; NaN where a row has no label). For each labeling, the
    column is split by class and the *same* classify-and-stat is applied to each
    subset, emitted under ``by_<labeling>`` keyed by class label. Rows whose
    label is missing are excluded from every class (but still counted in the
    overall stats above). The stat kind (numeric/categorical) is decided ONCE on
    the full column so every class uses the same shape.
    """
    overall = _classify_and_stat(series)
    if not groupings:
        return overall
    is_numeric = overall.get("kind") == "numeric"
    for labeling_name, labels in groupings.items():
        by_class: Dict[str, object] = {}
        # Align labels to this column's index; drop rows with no label.
        aligned = labels.reindex(series.index)
        valid_mask = aligned.notna()
        for class_label, idx in aligned[valid_mask].groupby(aligned[valid_mask]).groups.items():
            subset = series.loc[idx]
            if is_numeric:
                by_class[str(class_label)] = _numeric_stats(subset)
            else:
                by_class[str(class_label)] = _categorical_stats(subset)
        overall[f"by_{labeling_name}"] = by_class
    return overall


# Column names that identify the row's protein, in priority order. The first one
# present becomes the canonical ``ID`` (tools vary: ``ID`` / lowercase ``id`` from EE
# seq-only / ``fa_id`` from SoluProt). Lets the labeling join + per-design drop work
# regardless of which a tool emitted.
_ID_ALIASES = ("ID", "id", "Id", "fa_id", "reference_id")


def _normalize_id_column(df: "pd.DataFrame") -> "pd.DataFrame":
    """Rename whichever ID-alias column is present to the canonical ``ID``.

    Leaves df unchanged if ``ID`` already exists. Surplus alias columns (and a
    pure row-index ``runtime_id``) are dropped downstream via ``_DROP_COLUMNS``.
    """
    if "ID" in df.columns:
        return df
    for alias in _ID_ALIASES:
        if alias in df.columns:
            return df.rename(columns={alias: "ID"})
    return df


def aggregate_csv(
    csv_path: str,
    labelings: Optional[Dict[str, Dict[str, str]]] = None,
) -> Dict[str, object]:
    """Aggregate one metric CSV into ``{n_rows, source_csv, columns: {...}}``.

    ``labelings`` maps ``labeling_name -> {reference_id: label}``. When given,
    each numeric/categorical column additionally carries a ``by_<labeling>``
    block with the same stats computed per class (join: this CSV's ``ID`` vs the
    labeling's ``reference_id``).
    """
    df = pd.read_csv(csv_path)
    df = _normalize_id_column(df)
    # Build per-labeling label Series aligned to df's rows (by the ID column).
    groupings: Dict[str, pd.Series] = {}
    if labelings and "ID" in df.columns:
        ids = df["ID"].astype(str)
        for labeling_name, mapping in labelings.items():
            groupings[labeling_name] = ids.map(mapping)
    columns: Dict[str, object] = {}
    for col in df.columns:
        if col in _DROP_COLUMNS:
            continue
        columns[col] = _column_stats(df[col], groupings or None)
    return {
        "n_rows": int(len(df)),
        "source_csv": os.path.basename(csv_path),
        "columns": columns,
    }


def _metric_from_filename(base: str) -> Tuple[str, bool]:
    """Derive (metric_name, is_comparative) from a CSV basename.

    Prefers a canonical name from METRIC_SUFFIXES; flags comparative metrics; and
    for an unknown (newly-added) metric strips a known input prefix so it still
    bands without a code change here.
    """
    stem = base[:-len(".csv")] if base.endswith(".csv") else base
    for metric, suffix in sorted(METRIC_SUFFIXES.items(), key=lambda kv: len(kv[1]), reverse=True):
        if stem == suffix or stem.endswith("_" + suffix):
            return metric, False
    for suffix in sorted(_COMPARATIVE_SUFFIXES, key=len, reverse=True):
        if stem == suffix or stem.endswith("_" + suffix):
            return suffix, True
    for prefix in _KNOWN_PREFIXES:
        if stem.startswith(prefix):
            return stem[len(prefix):], False
    return stem, False


def discover_csvs(input_dir: str) -> Dict[str, str]:
    """Map metric-name -> CSV path for EVERY metric CSV present in input_dir.

    Generic: bands every ``*.csv`` keyed by ID, so a newly-added metric is picked up
    with no code change here. Skips labeling files (``*_labels.csv``) and the
    inherently comparative metrics (``_COMPARATIVE_SUFFIXES``). If several CSVs map to
    one metric, the lexicographically last is used and a note is printed.
    """
    found: Dict[str, str] = {}
    for path in sorted(glob.glob(os.path.join(input_dir, "*.csv"))):
        base = os.path.basename(path)
        if base.endswith("_labels.csv"):
            continue  # labeling files are consumed as by_<labeling> strata, not banded
        metric, is_comparative = _metric_from_filename(base)
        if is_comparative:
            print(f"[skip] {base}: '{metric}' is comparative (similarity to a reference set) -> no natural band")
            continue
        if metric in found:
            print(f"[warn] {metric}: multiple CSVs match, using last: {found[metric]} -> {path}")
        found[metric] = path
    return found


def load_label_file(path: str, labeling_name: Optional[str] = None) -> Dict[str, str]:
    """Load a ``reference_id,label`` (2-column) CSV into a ``{id: label}`` dict.

    Header-agnostic: the first column is the reference id, the second is the
    label, regardless of their header names (so ``reference_id,label`` and
    ``ID,domain_architecture`` both work). Blank labels are dropped. Duplicate
    ids keep the last occurrence.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    if df.shape[1] < 2:
        raise SystemExit(
            f"Label file {path!r} needs >=2 columns (reference_id,label); "
            f"found {list(df.columns)}."
        )
    id_col, label_col = df.columns[0], df.columns[1]
    mapping: Dict[str, str] = {}
    for ref_id, label in zip(df[id_col].astype(str), df[label_col].astype(str)):
        ref_id = ref_id.strip()
        label = label.strip()
        if ref_id and label != "":
            mapping[ref_id] = label
    return mapping


def labeling_name_from_path(path: str) -> str:
    """Derive a labeling name from a file path (basename minus extension)."""
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    # Strip a trailing ``_labels`` so first_cyclization_labels.csv -> first_cyclization.
    if stem.endswith("_labels"):
        stem = stem[: -len("_labels")]
    return stem


def build_reference_stats(
    input_dir: str,
    *,
    explicit: Optional[Dict[str, str]] = None,
    labelings: Optional[Dict[str, Dict[str, str]]] = None,
) -> Dict[str, object]:
    """Build the full reference-stats dict from a dir of per-metric CSVs.

    ``explicit`` optionally maps metric-name -> csv path, overriding/augmenting
    suffix-based discovery (used when a CSV lives elsewhere or is renamed).
    ``labelings`` maps labeling-name -> {reference_id: label}; when given, every
    metric column also gets a ``by_<labeling>`` per-class stratification block.
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
        out[metric] = aggregate_csv(path, labelings=labelings)
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
        help="Reference set label embedded in the JSON metadata (default: marts_db). "
        "Keep this version-free; record the version via --release_date / --structure_source.",
    )
    parser.add_argument(
        "--structure_source",
        default=None,
        help="Structure-prediction source backing the structure metrics (e.g. 'esmfold', "
        "'af3'). Recorded in the JSON metadata. Sequence metrics are source-independent.",
    )
    parser.add_argument(
        "--release_date",
        default=None,
        help="MARTS-DB release/version (e.g. '2026-06-12') this reference was built from. "
        "Recorded in the JSON metadata so the version travels WITH the data and the "
        "filename can stay version-free (stable across releases).",
    )
    parser.add_argument(
        "--group_by",
        action="append",
        default=[],
        metavar="LABEL_FILE[:NAME]",
        help="Add per-class stratification: a 'reference_id,label' (2-column) CSV "
        "label file. Each metric column then also carries a 'by_<NAME>' block "
        "with the same stats computed per class. NAME defaults to the file's "
        "basename (trailing '_labels' stripped). Repeatable for multiple labelings.",
    )
    parser.add_argument(
        "--group_by_column",
        action="append",
        default=[],
        metavar="METRIC:COLUMN[:NAME]",
        help="Add per-class stratification using a column of one of the metric "
        "CSVs in input_dir as the labeling (joined on ID). E.g. "
        "'domain_composition:domain_architecture'. NAME defaults to COLUMN. "
        "Repeatable.",
    )
    args = parser.parse_args()

    # Assemble labelings: external label files + columns lifted from metric CSVs.
    labelings: Dict[str, Dict[str, str]] = {}
    for spec in args.group_by:
        path, _, name = spec.partition(":")
        name = name or labeling_name_from_path(path)
        labelings[name] = load_label_file(path)
        print(f"[group_by] {name}: {len(labelings[name])} labels from {path}")
    if args.group_by_column:
        csvs = discover_csvs(args.input_dir)
        for spec in args.group_by_column:
            parts = spec.split(":")
            metric, column = parts[0], parts[1]
            name = parts[2] if len(parts) > 2 else column
            if metric not in csvs:
                raise SystemExit(
                    f"--group_by_column {spec!r}: metric {metric!r} not found in "
                    f"{args.input_dir!r}. Available: {sorted(csvs)}"
                )
            df = pd.read_csv(csvs[metric], dtype=str, keep_default_na=False)
            if "ID" not in df.columns or column not in df.columns:
                raise SystemExit(
                    f"--group_by_column {spec!r}: CSV {csvs[metric]!r} needs an "
                    f"'ID' column and a {column!r} column; has {list(df.columns)}."
                )
            mapping = {
                str(i).strip(): str(c).strip()
                for i, c in zip(df["ID"], df[column])
                if str(i).strip() and str(c).strip() != ""
            }
            labelings[name] = mapping
            print(f"[group_by] {name}: {len(mapping)} labels from {metric}:{column}")

    stats = build_reference_stats(args.input_dir, labelings=labelings or None)
    document = {
        "reference_set": args.reference_name,
        "structure_source": args.structure_source,
        "marts_db_release": args.release_date,
        "input_dir": os.path.abspath(args.input_dir),
        "labelings": sorted(labelings) if labelings else [],
        "metrics": _json_safe(stats),
    }
    with open(args.output, "w") as fh:
        json.dump(document, fh, indent=2, sort_keys=False)
        fh.write("\n")
    n_metrics = len(stats)
    print(f"\nWrote {n_metrics} metric(s) to {args.output}")


if __name__ == "__main__":
    main()
