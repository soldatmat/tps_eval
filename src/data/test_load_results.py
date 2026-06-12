from __future__ import annotations

"""Self-contained tests for load_results.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/data && python test_load_results.py
or under pytest:
    cd src/data && python -m pytest test_load_results.py -q

These tests use synthetic in-memory DataFrames; no CSVs / fasta files are read.
They lock in the empty-frames robustness fix (a structure-only run with no
sequence-branch outputs must not crash reduce() with an empty iterable) and that
the normal multi-frame outer-join happy path is unchanged.
"""

import sys
from pathlib import Path

# Put src/ on the path so the package-style imports inside load_results
# (`from data.sequences import ...`) resolve regardless of the cwd.
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pandas as pd

from data.load_results import _outer_join_on_id
from plot.boxplot_comparison import boxplot_comparison
from plot.categorical_comparison import categorical_comparison
from plot.density_comparison import density_comparison


def test_outer_join_empty_returns_id_frame():
    """Empty frames list -> empty DataFrame carrying the 'ID' key, no crash."""
    out = _outer_join_on_id([])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["ID"]
    assert len(out) == 0


def test_outer_join_single_frame_passthrough():
    df = pd.DataFrame({"ID": ["a", "b"], "x": [1.0, 2.0]})
    out = _outer_join_on_id([df])
    assert out is df  # single-frame fast path returns the frame untouched


def test_outer_join_happy_path_merges():
    """The normal non-empty case still outer-merges on ID (unchanged behavior)."""
    left = pd.DataFrame({"ID": ["a", "b"], "x": [1.0, 2.0]})
    right = pd.DataFrame({"ID": ["b", "c"], "y": [3.0, 4.0]})
    out = _outer_join_on_id([left, right])
    assert set(out.columns) == {"ID", "x", "y"}
    assert set(out["ID"]) == {"a", "b", "c"}
    # b has both x and y; a only x; c only y.
    by_id = out.set_index("ID")
    assert by_id.loc["b", "x"] == 2.0 and by_id.loc["b", "y"] == 3.0
    assert pd.isna(by_id.loc["a", "y"])
    assert pd.isna(by_id.loc["c", "x"])


def test_boxplot_skips_empty_target(capsys=None):
    """boxplot_comparison must skip (not KeyError) when no df has the target."""
    empty = pd.DataFrame(columns=["ID"])
    # Must not raise; produces no file (save_dir omitted).
    boxplot_comparison(
        ["__x__.fasta"], ["gen"], ["dodgerblue"], "soluprot", dfs=[empty]
    )


def test_density_skips_empty_target():
    """density_comparison must skip (not KeyError) when no df has the target."""
    empty = pd.DataFrame(columns=["ID"])
    density_comparison(
        ["__x__.fasta"], ["gen"], ["dodgerblue"], "soluprot", dfs=[empty]
    )


def test_categorical_skips_empty_target():
    """categorical_comparison must skip gracefully when no values found."""
    empty = pd.DataFrame(columns=["ID"])
    categorical_comparison(
        ["__x__.fasta"], ["gen"], ["dodgerblue"], "some_bool", dfs=[empty]
    )


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
