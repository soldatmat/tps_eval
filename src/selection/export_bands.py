"""export_bands — turn a reference-stats JSON into a band_filter ``bands_file``.

The reference-stats pipeline (scripts/compute_reference_stats.sh ->
src/reference_stats/aggregate_reference_stats.py) emits, for every metric column, a
percentile summary of the natural-TPS distribution — overall AND (with --group_by) a
``by_<labeling>`` block per class. This helper RESOLVES those percentiles into the
concrete ``[min, max]`` bands that band_filter consumes, so a funnel can pull its
reference bands straight from the natural distribution instead of hard-coding them.

Per-architecture bands (the capability the production funnel needed) = run the aggregator
with a single-vs-two-domain labeling (``--group_by domain.csv`` mapping reference_id ->
single|two), then export with ``--by domain`` here: each metric's band becomes a
``{"by": <labeling>, "<cat>": {min,max}}`` block keyed by the design's architecture.

Band edges are configurable per metric via a spec ({metric: {lo: <pctl|null>, hi: <pctl|null>}},
e.g. ``{"lo": "p25", "hi": "p75"}`` for the IQR, or ``{"lo": "p5", "hi": null}`` for a
one-sided floor). A default (lo/hi percentile pair) applies to any metric not named.
Output is band_filter's format: ``{"metrics": {metric: LEAF}}``.
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, Optional


def _edge(stats: dict, pctl: Optional[str]):
    if pctl is None:
        return None
    if pctl not in stats:
        raise ValueError(f"percentile '{pctl}' not in stats (have {sorted(stats)}).")
    return stats[pctl]


def _leaf_from_stats(stats: dict, lo: Optional[str], hi: Optional[str]) -> dict:
    leaf = {}
    lo_v, hi_v = _edge(stats, lo), _edge(stats, hi)
    if lo_v is not None:
        leaf["min"] = lo_v
    if hi_v is not None:
        leaf["max"] = hi_v
    return leaf


def export_bands(ref_stats: dict, metrics_spec: Dict[str, dict],
                 default_edges: dict, by: Optional[str] = None) -> dict:
    """Build a band_filter ``{"metrics": {...}}`` dict from a reference-stats JSON.

    ``metrics_spec`` maps metric -> {"lo": <pctl|null>, "hi": <pctl|null>}; ``default_edges``
    is the fallback pair. If ``by`` is given, each metric emits a per-category band block
    from the reference-stats ``by_<by>`` stratum (categories with no stratum are skipped).
    """
    ref_metrics = ref_stats.get("metrics", {})
    # Flatten reference-stats: metric-column -> its stats dict (+ its by_<...> strata).
    col_stats: Dict[str, dict] = {}
    for _tool, block in ref_metrics.items():
        for col, cstats in block.get("columns", {}).items():
            col_stats[col] = cstats

    out_metrics: Dict[str, dict] = {}
    for metric, edges in {**{m: default_edges for m in metrics_spec}, **metrics_spec}.items():
        if metric not in col_stats:
            print(f"  [export_bands] '{metric}' not in reference stats — skipped.", file=sys.stderr)
            continue
        lo, hi = edges.get("lo", default_edges.get("lo")), edges.get("hi", default_edges.get("hi"))
        cstats = col_stats[metric]
        if by:
            strata = cstats.get(f"by_{by}", {})
            leaf = {"by": by}
            for cat, cat_stats in strata.items():
                if cat_stats.get("kind") == "numeric":
                    leaf[cat] = _leaf_from_stats(cat_stats, lo, hi)
            if len(leaf) > 1:
                out_metrics[metric] = leaf
        else:
            out_metrics[metric] = _leaf_from_stats(cstats, lo, hi)
    return {"metrics": out_metrics, "_source": ref_stats.get("structure_source"),
            "_reference_set": ref_stats.get("reference_set")}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref_stats", required=True, help="Reference-stats JSON (aggregate output).")
    p.add_argument("--output", required=True, help="Output band_filter bands_file JSON.")
    p.add_argument("--metrics", nargs="+", required=True,
                   help="Metrics to band. Each 'metric' uses the default edges, or "
                        "'metric:lo,hi' to override (e.g. catalytic_pocket_volume:p25,p75 or "
                        "pocket_depth:p25,none).")
    p.add_argument("--default_lo", default="p25", help="Default lower-edge percentile (or 'none').")
    p.add_argument("--default_hi", default="p75", help="Default upper-edge percentile (or 'none').")
    p.add_argument("--by", default=None,
                   help="Emit per-category bands from the reference-stats by_<BY> strata "
                        "(e.g. domain_architecture / a single-vs-two labeling).")
    args = p.parse_args()

    def _norm(v):
        return None if v is None or str(v).lower() in ("none", "null") else v

    default_edges = {"lo": _norm(args.default_lo), "hi": _norm(args.default_hi)}
    metrics_spec = {}
    for m in args.metrics:
        if ":" in m:
            name, edges = m.split(":", 1)
            lo, hi = (edges.split(",") + ["", ""])[:2]
            metrics_spec[name] = {"lo": _norm(lo), "hi": _norm(hi)}
        else:
            metrics_spec[m] = dict(default_edges)

    with open(args.ref_stats) as fh:
        ref = json.load(fh)
    bands = export_bands(ref, metrics_spec, default_edges, by=args.by)
    with open(args.output, "w") as fh:
        json.dump(bands, fh, indent=2)
    print(f"[export_bands] wrote {len(bands['metrics'])} metric band(s) -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
