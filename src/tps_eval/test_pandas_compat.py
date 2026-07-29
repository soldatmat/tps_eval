"""Unit tests for the pandas major-version shims.

Regression: `groupby.idxmax()` on an ALL-NA group returned NaN under pandas 2 but
RAISES under pandas 3. The foldseek best-hit reducers (structure_alignment /
domain_alignment -> structural_identity, the TM-score fold gate) depend on the NaN,
so on pandas 3 the whole reduction died on a single unscorable query.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tps_eval.pandas_compat import group_idxmax_skipna


def _frame():
    return pd.DataFrame(
        {
            "query": ["a", "a", "b", "b", "c", "c"],
            "lddt": [0.1, 0.9, np.nan, 0.5, np.nan, np.nan],  # c is ALL-NA
        }
    )


def test_picks_the_max_row_per_group():
    df = _frame()
    out = group_idxmax_skipna(df.groupby("query"), "lddt")
    assert out.loc["a"] == 1  # 0.9 beats 0.1
    assert out.loc["b"] == 3  # the only non-NA in b
    assert df.loc[out.loc["a"], "lddt"] == 0.9


def test_all_na_group_yields_nan_instead_of_raising():
    df = _frame()
    out = group_idxmax_skipna(df.groupby("query"), "lddt")
    assert pd.isna(out.loc["c"])
    assert list(out.index) == ["a", "b", "c"]  # every group present


def test_every_group_all_na():
    df = pd.DataFrame({"query": ["a", "a", "b"], "lddt": [np.nan] * 3})
    out = group_idxmax_skipna(df.groupby("query"), "lddt")
    assert list(out.index) == ["a", "b"]
    assert out.isna().all()


def test_matches_native_idxmax_when_no_group_is_all_na():
    """Where plain idxmax is legal, the shim must agree with it exactly."""
    df = pd.DataFrame(
        {"query": ["a", "a", "b", "b"], "lddt": [0.2, np.nan, 0.7, 0.9]}
    )
    grouped = df.groupby("query")
    native = grouped["lddt"].idxmax()
    shimmed = group_idxmax_skipna(grouped, "lddt")
    assert list(shimmed.index) == list(native.index)
    assert list(shimmed) == list(native)
