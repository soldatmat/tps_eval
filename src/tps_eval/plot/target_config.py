from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np


def resolve_range(
    target: str,
    min_val: Dict[str, float],
    max_val: Dict[str, float],
    all_data: Sequence[Sequence[float]],
) -> Tuple[float, float]:
    """Axis (min, max) for a target.

    Uses the fixed bounds from the `MIN_VAL`/`MAX_VAL` constants when defined;
    otherwise derives a padded range from the pooled data so newly-added metrics
    (which have no hardcoded scale) plot sensibly. Falls back to (0, 1) when the
    target has no finite data at all.
    """
    if target in min_val and target in max_val:
        return min_val[target], max_val[target]

    pooled = np.asarray(
        [v for series in all_data for v in series if np.isfinite(v)],
        dtype=float,
    )
    if pooled.size == 0:
        return 0.0, 1.0

    lo = float(pooled.min())
    hi = float(pooled.max())
    if lo == hi:
        # Degenerate single-value distribution: pad symmetrically.
        pad = abs(lo) * 0.05 if lo != 0.0 else 0.5
        return lo - pad, hi + pad
    pad = (hi - lo) * 0.05
    return lo - pad, hi + pad


def auto_ticks(min_val: float, max_val: float, n: int = 11) -> np.ndarray:
    """Evenly-spaced tick array spanning [min_val, max_val] for auto-ranged
    targets (those without an explicit `TICKS` entry)."""
    return np.round(np.linspace(min_val, max_val, n), 6)
