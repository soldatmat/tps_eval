from __future__ import annotations

"""Self-contained tests for band_filter.py (keep designs within reference bands).

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_band_filter.py
or under pytest:
    cd src/selection && python -m pytest test_band_filter.py -q

Synthetic in-memory DataFrames + a tiny temp bands_file JSON. Checks the inclusive
[min,max] boundaries, one-sided bands, the per-architecture (`by`) block incl. the
skip-uncovered-category rule, the missing-value-fails rule, bands_file loading with
inline override, and the unknown-column guards.
"""

import json
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from band_filter import apply_band_filter  # noqa: E402


def _ids(df):
    return set(df["ID"])


def _df():
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "v": [10.0, 20.0, 30.0, np.nan],
        "arch": ["single", "single", "two", "two"],
        "vol": [800.0, 200.0, 1000.0, 1500.0],
    })


def test_inclusive_boundaries():
    out, rep = apply_band_filter(_df(), {"v": {"min": 10.0, "max": 30.0}})
    assert _ids(out) == {"a", "b", "c"}  # both endpoints inside; d (NaN) fails
    assert rep["n_in"] == 4 and rep["n_pass"] == 3


def test_one_sided_min_only():
    out, _ = apply_band_filter(_df(), {"v": {"min": 20.0}})
    assert _ids(out) == {"b", "c"}


def test_one_sided_max_only():
    out, _ = apply_band_filter(_df(), {"v": {"max": 20.0}})
    assert _ids(out) == {"a", "b"}


def test_missing_value_fails():
    out, _ = apply_band_filter(_df(), {"v": {"min": 0.0}})
    assert "d" not in _ids(out)  # NaN cannot be shown in-range


def test_per_architecture_bands():
    out, _ = apply_band_filter(_df(), {
        "vol": {"by": "arch",
                "single": {"min": 617, "max": 1377},
                "two": {"min": 326, "max": 1016}},
    })
    # single: a(800) in [617,1377] pass, b(200) fail; two: c(1000) in [326,1016] pass,
    # d(1500) fail.
    assert _ids(out) == {"a", "c"}


def test_per_architecture_uncovered_category_passes():
    # Provide a band only for 'single'; 'two' rows have no band -> not filtered (pass).
    out, _ = apply_band_filter(_df(), {
        "vol": {"by": "arch", "single": {"min": 617, "max": 1377}},
    })
    # single: a passes, b fails; two: c,d both auto-pass (no band).
    assert _ids(out) == {"a", "c", "d"}


def test_bands_file_with_inline_override():
    with tempfile.TemporaryDirectory() as d:
        bf = os.path.join(d, "bands.json")
        with open(bf, "w") as fh:
            json.dump({"metrics": {"v": {"min": 0.0, "max": 100.0}}}, fh)
        # Inline metrics override the file's band for the same metric.
        out, _ = apply_band_filter(_df(), {"v": {"min": 25.0}}, bands_file=bf)
    assert _ids(out) == {"c"}  # inline min=25 wins over file's [0,100]


def test_bands_file_only():
    with tempfile.TemporaryDirectory() as d:
        bf = os.path.join(d, "bands.json")
        with open(bf, "w") as fh:
            json.dump({"metrics": {"vol": {"min": 500.0, "max": 1200.0}}}, fh)
        out, _ = apply_band_filter(_df(), {}, bands_file=bf)
    assert _ids(out) == {"a", "c"}  # 800 and 1000 in [500,1200]


def test_unknown_metric_raises():
    try:
        apply_band_filter(_df(), {"nope": {"min": 0}})
    except ValueError as e:
        assert "unknown column" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown metric")


def test_unknown_by_column_raises():
    try:
        apply_band_filter(_df(), {"vol": {"by": "missing_col", "single": {"min": 1}}})
    except ValueError as e:
        assert "by-column" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown by-column")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
