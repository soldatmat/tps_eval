"""band_filter — keep designs whose metrics fall within reference (natural-TPS) bands.

A band is an inclusive ``[min, max]`` per metric. The genuinely new capability here (the
reason selection needed its own filter rather than the dashboard's) is PER-ARCHITECTURE
bands: single-domain TPS pockets are legitimately larger/deeper/more enclosed than
two-domain ones, so a single mixed band unfairly penalises one architecture. A metric's
band may therefore be conditioned on a categorical column (e.g. ``domain_architecture`` ∈
{single, two}), each architecture getting its own ``[min, max]``.

Band spec (``metrics``): ``{metric: LEAF}`` where LEAF is either
  * ``{"min": x, "max": y}``  (either bound optional -> one-sided), or
  * ``{"by": "<col>", "<value>": {"min", "max"}, ...}``  (per-category bands).
Bands may also be loaded from a resolved-band JSON (``bands_file``) produced by the
reference-stats pipeline; inline ``metrics`` override/augment it.

A design with a missing metric value FAILS that metric's band (it cannot be shown to be in
range). A per-``by`` band whose category value is absent from the spec is SKIPPED for that
design (not applied), with a one-time warning.
"""
from __future__ import annotations

import json
from typing import Dict, Tuple

import pandas as pd


def _band_mask(vals: pd.Series, band: dict) -> pd.Series:
    """Vectorised in-range mask for a numeric series against one ``{min?,max?}`` band.
    A missing value is out-of-range (cannot be shown to be within)."""
    ok = vals.notna()
    if band.get("min") is not None:
        ok &= vals >= band["min"]
    if band.get("max") is not None:
        ok &= vals <= band["max"]
    return ok


def apply_band_filter(df: pd.DataFrame, metrics: Dict[str, dict],
                      bands_file: str = None) -> Tuple[pd.DataFrame, Dict]:
    """Keep rows within every metric's band; add ``band_pass`` and drop failers."""
    resolved: Dict[str, dict] = {}
    if bands_file:
        with open(bands_file) as fh:
            resolved.update(json.load(fh).get("metrics", {}))
    resolved.update(metrics or {})

    mask = pd.Series(True, index=df.index)
    per_metric = []
    for metric, leaf in resolved.items():
        if metric not in df.columns:
            raise ValueError(f"band_filter references unknown column '{metric}'.")
        vals = pd.to_numeric(df[metric], errors="coerce")
        if "by" in leaf:
            by_col = leaf["by"]
            if by_col not in df.columns:
                raise ValueError(f"band_filter '{metric}' by-column '{by_col}' not present.")
            # Rows whose category has no band are NOT filtered on this metric (pass).
            mmask = pd.Series(True, index=df.index)
            cats = df[by_col].astype("string")
            covered = {k for k in leaf if k != "by"}
            for cat in cats.dropna().unique():
                if cat not in covered:
                    print(f"  [band_filter] '{metric}': no band for {by_col}='{cat}' "
                          f"-> not applied to those rows.")
                    continue
                sel = cats == cat
                mmask.loc[sel] = _band_mask(vals[sel], leaf[cat])
        else:
            mmask = _band_mask(vals, leaf)
        per_metric.append({"metric": metric, "passed": int(mmask.sum())})
        mask &= mmask

    out = df.copy()
    out["band_pass"] = mask.values
    report = {"op": "band_filter", "n_in": len(df), "n_pass": int(mask.sum()),
              "metrics": per_metric}
    out = out[out["band_pass"]].drop(columns=["band_pass"])
    return out, report
