from __future__ import annotations

"""Self-contained tests for score.py (weighted z-sum ranking within a group).

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_score.py
or under pytest:
    cd src/selection && python -m pytest test_score.py -q

Synthetic in-memory DataFrames whose z-scores are known in closed form, so the
score arithmetic is checked exactly. Covers within-group z (different group scales
compared fairly), direction sign-flip, the constant-column divide-by-zero guard
(contributes 0, not NaN), the missing-value -> NaN-score / last-rank contract, the
global (no-group) pool, and the unknown-column guards.
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from score import _zscore, apply_score  # noqa: E402


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b}"


def test_zscore_closed_form():
    z = _zscore(pd.Series([1.0, 2.0, 3.0]))
    # mean 2, population std (ddof=0) = sqrt(2/3).
    std = np.sqrt(2.0 / 3.0)
    _approx(z.iloc[0], (1 - 2) / std)
    _approx(z.iloc[2], (3 - 2) / std)


def test_direction_flip():
    df = pd.DataFrame({"ID": ["a", "b"], "nll": [1.0, 3.0]})
    higher, _ = apply_score(df, [{"col": "nll", "weight": 1, "direction": "higher"}])
    lower, _ = apply_score(df, [{"col": "nll", "weight": 1, "direction": "lower"}])
    # higher: a<b; lower flips it: a>b.
    assert higher.set_index("ID").loc["a", "score"] < higher.set_index("ID").loc["b", "score"]
    assert lower.set_index("ID").loc["a", "score"] > lower.set_index("ID").loc["b", "score"]


def test_within_group_scaling():
    # Two groups on very different scales; within-group z puts them on equal footing.
    df = pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "g": ["big", "big", "small", "small"],
        "v": [1000.0, 2000.0, 1.0, 2.0],
    })
    out, rep = apply_score(df, [{"col": "v", "weight": 1, "direction": "higher"}],
                           zscore_within="g")
    s = out.set_index("ID")["score"]
    # Top of each group has the same +z; bottom the same -z.
    _approx(s["a"], s["c"])
    _approx(s["b"], s["d"])
    assert rep["zscore_within"] == "g"


def test_constant_column_contributes_zero_not_nan():
    df = pd.DataFrame({"ID": ["a", "b"], "c": [5.0, 5.0], "v": [1.0, 2.0]})
    out, _ = apply_score(df, [
        {"col": "c", "weight": 1, "direction": "higher"},   # constant -> 0 contribution
        {"col": "v", "weight": 1, "direction": "higher"},
    ])
    s = out.set_index("ID")["score"]
    assert s.notna().all()  # constant column must not wipe out the score
    assert s["b"] > s["a"]


def test_missing_value_gives_nan_score_and_last_rank():
    df = pd.DataFrame({"ID": ["a", "b", "c"], "v": [1.0, np.nan, 3.0]})
    out, rep = apply_score(df, [{"col": "v", "weight": 1, "direction": "higher"}])
    s = out.set_index("ID")
    assert np.isnan(s.loc["b", "score"])
    assert np.isnan(s.loc["b", "score_rank"])  # NaN scores get no rank
    assert s.loc["c", "score_rank"] == 1.0     # best score ranks first
    assert rep["n_scored"] == 2


def test_weight_applied():
    df = pd.DataFrame({"ID": ["a", "b"], "v": [1.0, 3.0]})
    w1, _ = apply_score(df, [{"col": "v", "weight": 1, "direction": "higher"}])
    w2, _ = apply_score(df, [{"col": "v", "weight": 2, "direction": "higher"}])
    _approx(w2.set_index("ID").loc["b", "score"],
            2 * w1.set_index("ID").loc["b", "score"])


def test_global_pool_when_no_group():
    df = pd.DataFrame({"ID": ["a", "b", "c"], "v": [3.0, 2.0, 1.0]})
    out, _ = apply_score(df, [{"col": "v", "weight": 1, "direction": "higher"}])
    assert list(out.set_index("ID")["score_rank"].loc[["a", "b", "c"]]) == [1.0, 2.0, 3.0]


def test_unknown_group_column_raises():
    df = pd.DataFrame({"ID": ["a"], "v": [1.0]})
    try:
        apply_score(df, [{"col": "v"}], zscore_within="nope")
    except ValueError as e:
        assert "zscore_within" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown group column")


def test_unknown_term_column_raises():
    df = pd.DataFrame({"ID": ["a"], "v": [1.0]})
    try:
        apply_score(df, [{"col": "nope"}])
    except ValueError as e:
        assert "unknown column" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown term column")


def test_missing_group_key_retained_with_nan_score():
    # A row whose group key is NaN can't be z-scored within a class: it is RETAINED
    # (not silently dropped) with score=NaN / no rank, and counted via n_missing_group_key.
    df = pd.DataFrame({
        "ID": ["a", "b", "c"],
        "g": ["x", "x", np.nan],
        "v": [1.0, 2.0, 5.0],
    })
    out, rep = apply_score(df, [{"col": "v", "weight": 1, "direction": "higher"}],
                           zscore_within="g")
    s = out.set_index("ID")
    assert set(out["ID"]) == {"a", "b", "c"}      # NaN-group row retained, not dropped
    assert np.isnan(s.loc["c", "score"])
    assert np.isnan(s.loc["c", "score_rank"])
    assert rep["n_missing_group_key"] == 1
    assert s.loc["b", "score"] > s.loc["a", "score"]   # the valid group still scores normally


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
