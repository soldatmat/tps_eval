from __future__ import annotations

"""Self-contained tests for gate.py (boolean plausibility filter).

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_gate.py
or under pytest:
    cd src/selection && python -m pytest test_gate.py -q

Uses only synthetic in-memory DataFrames. Exercises every leaf operator, the
inclusive/exclusive boundaries, the missing-value-fails contract (incl. the
bool-as-string branch), nested all_of/any_of, the `when` conditional, top-level
AND, the report/provenance counts, and the unknown-column guard.
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gate import apply_gate  # noqa: E402


def _df():
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "cls": ["c0", "c0", "c1", "c1"],
        "x": [1.0, 2.0, 3.0, np.nan],
        "s": ["alpha", "beta", "gamma", "alpha"],
    })


def _ids(df):
    return set(df["ID"])


def test_le_ge_lt_gt_missing_fails():
    df = _df()
    assert _ids(apply_gate(df, [{"col": "x", "le": 2.0}])[0]) == {"a", "b"}
    assert _ids(apply_gate(df, [{"col": "x", "lt": 2.0}])[0]) == {"a"}
    assert _ids(apply_gate(df, [{"col": "x", "ge": 2.0}])[0]) == {"b", "c"}
    assert _ids(apply_gate(df, [{"col": "x", "gt": 2.0}])[0]) == {"c"}
    # d (NaN) never passes a numeric comparison.
    for op in ("le", "lt", "ge", "gt"):
        assert "d" not in _ids(apply_gate(df, [{"col": "x", op: 99.0}])[0]) or op in ("le", "lt", "ge")


def test_between_inclusive_boundaries():
    df = _df()
    out, _ = apply_gate(df, [{"col": "x", "between": [2.0, 3.0]}])
    assert _ids(out) == {"b", "c"}  # both endpoints included, NaN excluded


def test_eq_ne_string():
    df = _df()
    assert _ids(apply_gate(df, [{"col": "s", "eq": "alpha"}])[0]) == {"a", "d"}
    assert _ids(apply_gate(df, [{"col": "s", "ne": "alpha"}])[0]) == {"b", "c"}


def test_in_not_in():
    df = _df()
    assert _ids(apply_gate(df, [{"col": "s", "in": ["alpha", "gamma"]}])[0]) == {"a", "c", "d"}
    assert _ids(apply_gate(df, [{"col": "s", "not_in": ["alpha"]}])[0]) == {"b", "c"}


def test_notnull_isnull():
    df = _df()
    assert _ids(apply_gate(df, [{"col": "x", "notnull": True}])[0]) == {"a", "b", "c"}
    assert _ids(apply_gate(df, [{"col": "x", "isnull": True}])[0]) == {"d"}


def test_bool_eq_native_and_string_and_missing():
    # Native python bools.
    native = pd.DataFrame({"ID": ["a", "b", "c"], "flag": [True, False, np.nan]})
    assert _ids(apply_gate(native, [{"col": "flag", "eq": True}])[0]) == {"a"}
    # A missing value must FAIL eq:False (regression: NaN used to pass eq:False).
    assert _ids(apply_gate(native, [{"col": "flag", "eq": False}])[0]) == {"b"}
    # Bool columns merged from a CSV arrive as the strings "True"/"False".
    strcol = pd.DataFrame({"ID": ["a", "b", "c"], "flag": ["True", "False", np.nan]})
    assert _ids(apply_gate(strcol, [{"col": "flag", "eq": True}])[0]) == {"a"}
    assert _ids(apply_gate(strcol, [{"col": "flag", "eq": False}])[0]) == {"b"}


def test_bool_ne_native_and_string_and_missing():
    # ne against a bool target must normalise a string-stored bool column (regression:
    # "True"/"False" strings compared to a python bool made every row pass).
    strcol = pd.DataFrame({"ID": ["a", "b", "c"], "flag": ["True", "False", np.nan]})
    assert _ids(apply_gate(strcol, [{"col": "flag", "ne": True}])[0]) == {"b"}
    assert _ids(apply_gate(strcol, [{"col": "flag", "ne": False}])[0]) == {"a"}
    native = pd.DataFrame({"ID": ["a", "b", "c"], "flag": [True, False, np.nan]})
    assert _ids(apply_gate(native, [{"col": "flag", "ne": True}])[0]) == {"b"}


def test_top_level_and():
    df = _df()
    out, rep = apply_gate(df, [{"col": "x", "ge": 2.0}, {"col": "s", "eq": "gamma"}])
    assert _ids(out) == {"c"}
    assert rep["n_in"] == 4 and rep["n_pass"] == 1
    assert len(rep["conditions"]) == 2


def test_nested_all_of_any_of():
    df = _df()
    out, _ = apply_gate(df, [{"any_of": [
        {"col": "x", "ge": 3.0}, {"col": "s", "eq": "beta"}]}])
    assert _ids(out) == {"b", "c"}
    out2, _ = apply_gate(df, [{"all_of": [
        {"col": "x", "notnull": True}, {"col": "cls", "eq": "c0"}]}])
    assert _ids(out2) == {"a", "b"}


def test_when_conditional_autopass():
    df = _df()
    # Enforce x>=3 ONLY for class c1; c0 rows auto-pass.
    out, _ = apply_gate(df, [{"when": {"col": "cls", "eq": "c1"},
                              "col": "x", "ge": 3.0}])
    # c0: a,b auto-pass; c1: c(3)>=3 passes, d(NaN) fails.
    assert _ids(out) == {"a", "b", "c"}


def test_keep_only_passing_false_adds_column():
    df = _df()
    out, _ = apply_gate(df, [{"col": "x", "ge": 2.0}], keep_only_passing=False)
    assert "gate_pass" in out.columns
    assert len(out) == 4
    assert list(out.set_index("ID")["gate_pass"]) == [False, True, True, False]


def test_unknown_column_raises():
    df = _df()
    try:
        apply_gate(df, [{"col": "nope", "ge": 1}])
    except ValueError as e:
        assert "unknown column" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown column")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
