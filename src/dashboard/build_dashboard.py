"""Build a self-contained interactive HTML dashboard for the MARTS-DB reference
bands, with optional overlay of a generated-design batch.

The reference bands are the committed JSONs in ``src/reference_stats/``
(``marts_db_<structure_source>_metric_stats.json``). Each is a nested
``structure_source -> metric -> column -> stats`` tree where numeric columns
carry percentile bands (p1/p5/p25/median/p75/p95/p99 + min/max/mean/std) and a
``by_<labeling>`` breakdown, and categorical columns carry ``frequencies``.

A generated-design batch is the tps_eval pipeline's per-tool output: a directory
(or glob, or list) of ``<input>_<tool>.csv`` files (or a single merged CSV) keyed
by ``ID`` whose column names match the band column names exactly. Its raw
per-design values are overlaid on the same axis as the natural band.

Robust to missing bands: a design column with NO reference band (e.g. the
sequence / comparative metrics that have no natural band by design, or a
freshly-added metric) is still shown — as a design-only track, grouped under its
tool, so the batch can always be inspected even when some or all bands are absent.

Output is one self-contained HTML file (all CSS/JS/data inlined, no external
requests) so it is portable and viewable as a local file or a published artifact.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from typing import Dict, List, Optional, Tuple

try:
    from metric_info import METRIC_INFO, METRIC_CATEGORY, CATEGORY_ORDER
except ImportError:  # when imported as a module rather than run as a script
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from metric_info import METRIC_INFO, METRIC_CATEGORY, CATEGORY_ORDER

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REFERENCE_STATS_DIR = os.path.join(REPO_ROOT, "src", "reference_stats")
TEMPLATE_PATH = os.path.join(HERE, "template.html")
DATA_TOKEN = "/*__DASHBOARD_DATA__*/"

# Friendly labels for the structure sources.
SOURCE_LABELS = {
    "esmfold": "ESMFold",
    "af3": "AlphaFold3",
    "boltz2": "Boltz-2 (holo)",
    "boltz2_holo": "Boltz-2 (holo)",
}

# Map the JSON's `by_<labeling>` keys to short labeling names used in the UI.
LABELING_KEYS = {
    "by_substrate": "substrate",
    "by_first_cyclization": "first_cyclization",
    "by_domain_architecture": "domain_architecture",
}

# Numeric summary fields carried through to the UI (drop the rest to keep size down).
NUMERIC_FIELDS = (
    "count", "n_missing", "mean", "std", "min",
    "p1", "p5", "p25", "median", "p75", "p95", "p99", "max",
)

# Display order of metrics: fold quality first, then geometry / pocket / motif /
# composition / chemistry, then holo-only. Unknown metrics are appended.
METRIC_ORDER = [
    "plddt",
    "global_confidence",
    "interdomain_pae",
    "esm_pseudo_perplexity",
    "proteinmpnn_score",
    "radius_of_gyration",
    "active_site_geometry",
    "pocket_descriptors",
    "aromatic_lining",
    "diphosphate_sensor",
    "motif_pair_distance",
    "motif_structural_distance",
    "domain_composition",
    "aggregation",
    "ion_site_check",
    "substrate_positioning",
]

# How a design CSV filename's `_<tool>.csv` suffix maps to a metric-card label, for
# grouping design-only columns (those with no reference band). Specific multi-token
# suffixes that don't read well stripped are mapped explicitly; the rest fall back to
# the suffix token itself.
SUFFIX_METRIC = {
    "_motifs.csv": "motif_search",
    "_embedding_esm1b_min_embedding_distance.csv": "min_embedding_distance",
    "_embedding_esm1b_min_embedding_distance_self.csv": "min_embedding_distance_self",
}
_EXTRA_SUFFIXES = [
    "_motifs.csv", "_motif_pair_distance.csv",
    "_max_sequence_identity.csv", "_max_sequence_identity_self.csv",
    "_embedding_esm1b_min_embedding_distance.csv",
    "_embedding_esm1b_min_embedding_distance_self.csv",
    "_local_sequence_search.csv", "_local_sequence_search_self.csv",
    "_soluprot.csv", "_enzyme_explorer_sequence_only.csv", "_swissprot_search.csv",
    "_knn_label_transfer.csv", "_substrate_class.csv", "_sdr_divergence.csv",
    "_structural_identity.csv", "_foldseek_swissprot_search.csv",
    "_domain_structural_identity.csv", "_self_consistency.csv",
]
# Longest-first so e.g. `_max_sequence_identity_self.csv` wins over `_max_sequence_identity.csv`.
KNOWN_SUFFIXES = sorted(
    set([f"_{m}.csv" for m in METRIC_ORDER] + _EXTRA_SUFFIXES),
    key=len, reverse=True,
)

_MISSING = ("", "nan", "NA", "NaN", "None", "null")

# Identifier / per-design annotation columns that are never a banded metric (mirrors the
# aggregator's _DROP_COLUMNS). Dropped from design CSVs before overlay.
_DROP_DESIGN_COLUMNS = {
    "ID", "id", "Id", "fa_id", "reference_id", "runtime_id", "sequence",
    "ddxxd_motif", "ddxxd_start", "nse_dte_motif", "nse_dte_start", "best_template",
}


def _compact_numeric(col: dict) -> dict:
    return {k: col[k] for k in NUMERIC_FIELDS if k in col}


def _compact_categorical(col: dict) -> dict:
    out = {k: col[k] for k in ("count", "n_missing", "n_unique") if k in col}
    out["frequencies"] = col.get("frequencies", {})
    return out


def _compact_column(col: dict) -> dict:
    kind = col.get("kind", "numeric")
    out = _compact_categorical(col) if kind == "categorical" else _compact_numeric(col)
    out["kind"] = kind

    by = {}
    for json_key, short in LABELING_KEYS.items():
        if json_key in col and isinstance(col[json_key], dict):
            strata = {}
            for stratum, sub in col[json_key].items():
                if not isinstance(sub, dict):
                    continue
                strata[str(stratum)] = (
                    _compact_categorical(sub) if kind == "categorical"
                    else _compact_numeric(sub)
                )
            if strata:
                by[short] = strata
    if by:
        out["by"] = by
    return out


def _ordered_metrics(metric_names: List[str]) -> List[str]:
    known = [m for m in METRIC_ORDER if m in metric_names]
    extra = sorted(m for m in metric_names if m not in METRIC_ORDER)
    return known + extra


def load_band_source(path: str) -> dict:
    with open(path) as f:
        raw = json.load(f)
    metrics_in = raw.get("metrics", {})
    metrics_out = {}
    for mname in _ordered_metrics(list(metrics_in.keys())):
        m = metrics_in[mname]
        cols = m.get("columns", {})
        metrics_out[mname] = {
            "n_rows": m.get("n_rows"),
            "columns": {c: _compact_column(cols[c]) for c in cols},
        }
    src = raw.get("structure_source", "unknown")
    return {
        "structure_source": src,
        "label": SOURCE_LABELS.get(src, src),
        "n_rows": metrics_in[next(iter(metrics_in))].get("n_rows") if metrics_in else None,
        "release": raw.get("marts_db_release"),
        "reference_set": raw.get("reference_set"),
        "metrics": metrics_out,
    }


def _is_number(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def _parse_tool_label(path: str) -> str:
    base = os.path.basename(path)
    for suf in KNOWN_SUFFIXES:
        if base.endswith(suf):
            return SUFFIX_METRIC.get(suf, suf[1:-len(".csv")])
    return os.path.splitext(base)[0]


def _resolve_csv_paths(entries: List[str]) -> List[str]:
    """Expand a list of files / directories / globs into a deduped, ordered CSV list."""
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


def _band_column_kinds(sources: Dict[str, dict]) -> Dict[str, str]:
    """column name -> kind, unioned across all band sources."""
    kinds: Dict[str, str] = {}
    for src in sources.values():
        for m in src["metrics"].values():
            for cname, col in m["columns"].items():
                kinds.setdefault(cname, col.get("kind", "numeric"))
    return kinds


def load_design_batch(
    entries: List[str], sources: Dict[str, dict], name: Optional[str] = None
) -> Optional[dict]:
    """Load a generated-design batch as raw per-column value lists.

    Returns ``None`` if no usable CSVs are found. Columns are kept raw (floats for
    numeric, strings for categorical); ``col_tool`` records which tool each column
    came from so design-only columns (no reference band) can be grouped sensibly.
    """
    csv_paths = _resolve_csv_paths(entries)
    if not csv_paths:
        return None

    band_kind = _band_column_kinds(sources)
    merged: Dict[str, dict] = {}
    col_tool: Dict[str, str] = {}
    order: List[str] = []  # ID order of first appearance

    for p in csv_paths:
        tool = _parse_tool_label(p)
        with open(p, newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                continue
            # the row-id column varies by tool: ID / lowercase id / fa_id (SoluProt)
            id_col = next((a for a in ("ID", "id", "Id", "fa_id", "reference_id")
                           if a in reader.fieldnames), None)
            if id_col is None:
                continue
            data_cols = [c for c in reader.fieldnames if c not in _DROP_DESIGN_COLUMNS and c != id_col]
            for c in data_cols:
                col_tool.setdefault(c, tool)
            for r in reader:
                rid = r.get(id_col)
                if rid is None:
                    continue
                if rid not in merged:
                    merged[rid] = {}
                    order.append(rid)
                for c in data_cols:
                    if c not in merged[rid]:
                        merged[rid][c] = r.get(c)

    if not order:
        return None

    all_cols = list(col_tool.keys())
    # Decide each column's kind (band kind wins; else infer numeric-vs-categorical).
    values: Dict[str, List] = {}
    kinds: Dict[str, str] = {}
    for c in all_cols:
        raw = [merged[rid].get(c) for rid in order]
        nonnull = [v for v in raw if v is not None and v not in _MISSING]
        if c in band_kind:
            kind = band_kind[c]
        else:
            kind = "numeric" if nonnull and all(_is_number(v) for v in nonnull) else "categorical"
        kinds[c] = kind
        if kind == "numeric":
            values[c] = [float(v) if (v is not None and v not in _MISSING and _is_number(v)) else None
                         for v in raw]
        else:
            values[c] = [str(v) if (v is not None and v not in _MISSING) else None for v in raw]

    return {
        "name": name or (os.path.basename(os.path.dirname(os.path.abspath(csv_paths[0]))) or "designs"),
        "n": len(order),
        "synthetic": False,
        "ids": order,
        "values": values,
        "col_tool": col_tool,
        "col_kind": kinds,
        "n_files": len(csv_paths),
    }


def synthetic_design_batch(sources: Dict[str, dict], n: int = 28) -> dict:
    """A reproducible, clearly-labelled demo batch sampled around the ESMFold bands
    so the overlay mechanism is visible without a real pipeline run."""
    import random

    rng = random.Random(20260612)
    ref = sources.get("esmfold") or next(iter(sources.values()))
    headline = {
        "plddt": ["mean_plddt"],
        "global_confidence": ["ptm"],
        "radius_of_gyration": ["radius_of_gyration", "asphericity"],
        "active_site_geometry": ["carboxylate_convergence_radius", "catalytic_constellation_rmsd"],
        "pocket_descriptors": ["catalytic_pocket_volume", "pocket_hydrophobicity"],
        "aromatic_lining": ["n_pocket_aromatics", "aromatic_fraction"],
        "diphosphate_sensor": ["n_diphosphate_basic_residues"],
        "esm_pseudo_perplexity": ["esm_pseudo_perplexity"],
        "motif_structural_distance": ["motif_centroid_distance"],
        "domain_composition": ["n_domains"],
    }
    ids = [f"demo_design_{i + 1:03d}" for i in range(n)]
    values: Dict[str, List] = {}
    col_kind: Dict[str, str] = {}
    for mname, cols in headline.items():
        m = ref["metrics"].get(mname)
        if not m:
            continue
        for cname in cols:
            col = m["columns"].get(cname)
            if not col or col.get("kind") != "numeric":
                continue
            med, p25, p75 = col.get("median"), col.get("p25"), col.get("p75")
            if med is None or p25 is None or p75 is None:
                continue
            spread = max(p75 - p25, 1e-6)
            lo, hi = col.get("min"), col.get("max")
            out = []
            for _ in range(n):
                v = rng.gauss(med + 0.15 * spread, 0.9 * spread)
                if lo is not None:
                    v = max(v, lo - 0.3 * spread)
                if hi is not None:
                    v = min(v, hi + 0.3 * spread)
                out.append(round(v, 4))
            values[cname] = out
            col_kind[cname] = "numeric"
    return {
        "name": "synthetic demo batch", "n": n, "synthetic": True,
        "ids": ids, "values": values, "col_tool": {}, "col_kind": col_kind, "n_files": 0,
    }


def parse_design_specs(values: List[str]) -> List[Tuple[Optional[str], List[str]]]:
    """Parse ``--designs`` values into (name, [paths]) sets.

    Each value is one design SET: ``[name=]path[,path2,...]``. ``name=`` is optional
    (the part before the first comma is checked for ``=``); paths are comma-separated
    and each may be a file, directory, or glob.
    """
    specs: List[Tuple[Optional[str], List[str]]] = []
    for value in values:
        name: Optional[str] = None
        head = value.split(",", 1)[0]
        if "=" in head:
            name, _, rest = value.partition("=")
            name = name.strip() or None
            value = rest
        entries = [p.strip() for p in value.split(",") if p.strip()]
        if entries:
            specs.append((name, entries))
    return specs


def _design_only_groups(union_cols: Dict[str, Tuple[Optional[str], str]],
                        band_kind: Dict[str, str]) -> List[Tuple[str, dict]]:
    """Group design columns with NO reference band into metric cards.

    ``union_cols`` maps column -> (tool, kind) unioned across all design sets.
    Returns an ordered list of ``(metric_label, {column: {kind, band_missing}})``.
    """
    groups: Dict[str, dict] = {}
    for c, (tool, kind) in union_cols.items():
        if c in band_kind:
            continue  # band-backed -> rendered on its band card via the overlay
        groups.setdefault(tool or "design metrics", {})[c] = {
            "kind": kind, "band_missing": True,
        }
    ordered = [(m, groups[m]) for m in METRIC_ORDER if m in groups]
    ordered += [(m, groups[m]) for m in sorted(groups) if m not in METRIC_ORDER]
    return ordered


def build_data(source_paths: List[str],
               design_specs: Optional[List[Tuple[Optional[str], List[str]]]],
               demo: bool) -> dict:
    sources: Dict[str, dict] = {}
    for p in source_paths:
        src = load_band_source(p)
        sources[src["structure_source"]] = src

    release = next((s.get("release") for s in sources.values() if s.get("release")), None)
    reference_set = next(
        (s.get("reference_set") for s in sources.values() if s.get("reference_set")), "marts_db"
    )

    design_sets: List[dict] = []
    if design_specs:
        for name, entries in design_specs:
            ds = load_design_batch(entries, sources, name=name)
            if ds:
                design_sets.append(ds)
    elif demo and sources:
        design_sets.append(synthetic_design_batch(sources))

    # Append design-only metric cards (columns with no band) to every source so the
    # batch is always inspectable, even where bands are missing.
    if design_sets:
        band_kind = _band_column_kinds(sources)
        union_cols: Dict[str, Tuple[Optional[str], str]] = {}
        for ds in design_sets:
            for c in ds["values"]:
                if c not in union_cols:
                    union_cols[c] = (ds.get("col_tool", {}).get(c),
                                     ds.get("col_kind", {}).get(c, "numeric"))
        only_groups = _design_only_groups(union_cols, band_kind)
        if only_groups:
            if not sources:
                # No bands at all -> a pseudo "designs" source carrying every column.
                sources["designs"] = {
                    "structure_source": "designs", "label": "designs (no bands)",
                    "n_rows": sum(d["n"] for d in design_sets),
                    "release": release, "reference_set": reference_set,
                    "no_bands": True, "metrics": {},
                }
            for src in sources.values():
                for mlabel, cols in only_groups:
                    existing = src["metrics"].get(mlabel)
                    if existing is None:
                        src["metrics"][mlabel] = {"n_rows": None, "columns": dict(cols),
                                                  "design_only": True}
                    else:
                        for c, spec in cols.items():
                            existing["columns"].setdefault(c, spec)

    return {
        "reference_set": reference_set,
        "release": release,
        "labelings": ["substrate", "first_cyclization", "domain_architecture"],
        "sources": sources,
        "design_sets": design_sets,
        "metric_info": METRIC_INFO,
        "category_order": CATEGORY_ORDER,
        "metric_category": METRIC_CATEGORY,
    }


def render_html(data: dict) -> str:
    with open(TEMPLATE_PATH) as f:
        template = f.read()
    payload = "window.__DATA__ = " + json.dumps(data, allow_nan=False, separators=(",", ":")) + ";"
    if DATA_TOKEN not in template:
        raise RuntimeError(f"data injection token {DATA_TOKEN!r} not found in template")
    return template.replace(DATA_TOKEN, payload)


def _sanitize_for_json(obj):
    """Replace non-finite floats (NaN/inf) with None so the embedded JSON is strict."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize_for_json(v) for v in obj]
    return obj


def default_source_paths() -> List[str]:
    paths = []
    for src in ("esmfold", "af3", "boltz2"):
        p = os.path.join(REFERENCE_STATS_DIR, f"marts_db_{src}_metric_stats.json")
        if os.path.exists(p):
            paths.append(p)
    return paths


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--bands", nargs="*", default=None,
        help="Band JSON paths. Defaults to the committed MARTS-DB esmfold/af3/boltz2 JSONs. "
             "Pass an empty set / nonexistent paths to render a design-only view.",
    )
    ap.add_argument(
        "--designs", action="append", default=None, metavar="[NAME=]PATH[,PATH...]",
        help="A design SET to overlay. Repeatable — pass --designs once per set to overlay "
             "several sets (distinguished by dot-outline colour). Each value is "
             "'[name=]path[,path2,...]' where each path is a merged CSV, a directory, or a "
             "glob of the pipeline's *_<tool>.csv outputs (matched to bands by column name; "
             "columns with no band are still shown as design-only).",
    )
    ap.add_argument(
        "--demo", action="store_true",
        help="Overlay a synthetic demo batch (ignored if --designs is given).",
    )
    ap.add_argument(
        "--output", default=os.path.join(REPO_ROOT, "data", "dashboard", "marts_db_dashboard.html"),
        help="Output HTML path (default: data/dashboard/marts_db_dashboard.html).",
    )
    args = ap.parse_args(argv)

    source_paths = args.bands if args.bands is not None else default_source_paths()
    source_paths = [p for p in source_paths if os.path.exists(p)]

    design_specs = parse_design_specs(args.designs) if args.designs else None
    if not source_paths and not design_specs and not args.demo:
        ap.error("no band JSONs found and no --designs given; nothing to render.")

    data = build_data(source_paths, design_specs, demo=args.demo)
    if not data["sources"]:
        ap.error("nothing to render: no bands resolved and no design columns found.")
    data = _sanitize_for_json(data)
    html = render_html(data)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        f.write(html)

    n_sources = len(data["sources"])
    n_metrics = sum(len(s["metrics"]) for s in data["sources"].values())
    print(f"wrote {args.output}")
    print(f"  sources: {', '.join(data['sources']) or '(none)'}  ({n_sources} sources, {n_metrics} metric-tables)")
    for ds in data.get("design_sets", []):
        tag = " [synthetic demo]" if ds.get("synthetic") else f" from {ds.get('n_files', 0)} csv(s)"
        print(f"  design set '{ds['name']}': {ds['n']} designs{tag}")
    print(f"  size: {os.path.getsize(args.output) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
