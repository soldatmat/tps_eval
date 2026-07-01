"""score — rank designs by a weighted sum of z-scored metrics.

Only reliable, monotone, quality-tracking metrics belong here (fold confidence, designability,
low aggregation, interface confidence). Each term is z-scored WITHIN the grouping (so classes
with different metric scales are compared on equal footing), sign-flipped so "higher is better"
always, weighted, and summed into a ``score`` column. Designs with a missing value in ANY term
get ``score = NaN`` and rank last (they should have been gated out earlier).

Spec: ``terms = [{"col", "weight", "direction": "higher"|"lower"}, ...]`` and
``zscore_within`` = a grouping column name (or None for a single global pool).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _zscore(s: pd.Series) -> pd.Series:
    num = pd.to_numeric(s, errors="coerce")
    std = num.std(ddof=0)
    if not np.isfinite(std) or std == 0:
        # No spread -> the term carries no ranking signal; contribute 0 (not NaN, so a
        # constant column does not wipe out every design's score).
        return pd.Series(0.0, index=s.index).where(num.notna())
    return (num - num.mean()) / std


def apply_score(df: pd.DataFrame, terms: List[dict],
                zscore_within: Optional[str] = None,
                score_col: str = "score") -> Tuple[pd.DataFrame, Dict]:
    """Add ``score_col`` (+ ``<score_col>_rank`` within group) to ``df``."""
    out = df.copy()
    if zscore_within and zscore_within not in out.columns:
        raise ValueError(f"score zscore_within references unknown column '{zscore_within}'.")
    groups = out.groupby(zscore_within, sort=False) if zscore_within else [(None, out)]

    contrib = pd.Series(0.0, index=out.index)
    valid = pd.Series(True, index=out.index)
    for col_spec in terms:
        col = col_spec["col"]
        if col not in out.columns:
            raise ValueError(f"score term references unknown column '{col}'.")
        weight = float(col_spec.get("weight", 1.0))
        direction = col_spec.get("direction", "higher")
        sign = 1.0 if direction == "higher" else -1.0
        z = pd.Series(np.nan, index=out.index)
        for _, g in groups:
            z.loc[g.index] = _zscore(g[col])
        valid &= z.notna()
        contrib = contrib.add(sign * weight * z.fillna(0.0), fill_value=0.0)

    out[score_col] = contrib.where(valid)
    # Rank within group, best (highest score) = 1; NaN scores get no rank.
    if zscore_within:
        out[f"{score_col}_rank"] = (out.groupby(zscore_within, sort=False)[score_col]
                                    .rank(ascending=False, method="first"))
    else:
        out[f"{score_col}_rank"] = out[score_col].rank(ascending=False, method="first")
    report = {"op": "score", "n_in": len(df), "n_scored": int(valid.sum()),
              "terms": [{"col": t["col"], "weight": t.get("weight", 1.0),
                         "direction": t.get("direction", "higher")} for t in terms],
              "zscore_within": zscore_within}
    return out, report
