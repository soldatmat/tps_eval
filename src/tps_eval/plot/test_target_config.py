from __future__ import annotations

"""Self-contained tests for plot/target_config.py (axis-range / tick helpers).

Run from src/ so the package-style imports resolve like the plot runner does:
    cd src && python -m plot.test_target_config
or under pytest:
    cd src && python -m pytest plot/test_target_config.py -q

Pure numeric logic (no plotting): resolve_range prefers fixed bounds when both
are defined, else pads a data-derived range, handles the degenerate single-value
and all-empty/NaN cases; auto_ticks spans the window with n evenly-spaced ticks.
"""

import sys
from pathlib import Path


import numpy as np

from tps_eval.plot.target_config import auto_ticks, resolve_range


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_fixed_bounds_take_precedence():
    lo, hi = resolve_range("t", {"t": -0.01}, {"t": 1.01}, [[500.0, -500.0]])
    _approx(lo, -0.01)
    _approx(hi, 1.01)


def test_partial_bounds_fall_through_to_data():
    # Only MIN defined -> not "both defined" -> derive from data instead.
    lo, hi = resolve_range("t", {"t": 0.0}, {}, [[1.0, 2.0, 3.0]])
    pad = (3.0 - 1.0) * 0.05
    _approx(lo, 1.0 - pad)
    _approx(hi, 3.0 + pad)


def test_data_range_padded_5pct():
    lo, hi = resolve_range("t", {}, {}, [[0.0, 10.0], [5.0]])
    pad = (10.0 - 0.0) * 0.05
    _approx(lo, 0.0 - pad)
    _approx(hi, 10.0 + pad)


def test_nan_and_inf_ignored():
    lo, hi = resolve_range("t", {}, {}, [[float("nan"), 2.0, float("inf"), 4.0]])
    pad = (4.0 - 2.0) * 0.05
    _approx(lo, 2.0 - pad)
    _approx(hi, 4.0 + pad)


def test_all_empty_falls_back_to_unit():
    _approx_pair = resolve_range("t", {}, {}, [[], []])
    assert _approx_pair == (0.0, 1.0)


def test_all_nan_falls_back_to_unit():
    assert resolve_range("t", {}, {}, [[float("nan")], [float("inf")]]) == (0.0, 1.0)


def test_degenerate_single_value_nonzero():
    lo, hi = resolve_range("t", {}, {}, [[4.0, 4.0, 4.0]])
    pad = abs(4.0) * 0.05
    _approx(lo, 4.0 - pad)
    _approx(hi, 4.0 + pad)


def test_degenerate_single_value_zero():
    lo, hi = resolve_range("t", {}, {}, [[0.0, 0.0]])
    _approx(lo, -0.5)
    _approx(hi, 0.5)


def test_auto_ticks_span_and_count():
    ticks = auto_ticks(0.0, 1.0, n=11)
    assert ticks.shape == (11,)
    _approx(ticks[0], 0.0)
    _approx(ticks[-1], 1.0)
    assert np.all(np.diff(ticks) > 0)
    np.testing.assert_allclose(np.diff(ticks), 0.1, atol=1e-6)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
