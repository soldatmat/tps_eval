from __future__ import annotations

"""Self-contained tests for merge.py (per-tool CSVs -> one wide table keyed by ID).

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_merge.py
or under pytest:
    cd src/selection && python -m pytest test_merge.py -q

Writes tiny temp CSVs only. Locks in the merge conventions the dashboard shares:
per-cell FIRST-WINS over the UNION of IDs (a later file must not drop rows only it
has, and must not overwrite an earlier file's cell), id-alias -> canonical ID, the
`_self` suffixing, numeric coercion (all-or-nothing per column), the raw-feature-matrix
skip, and resolve_csv_paths ordering/dedup.
"""

import os
import sys
import tempfile

import numpy as np
import pandas as pd


from tps_eval.selection.merge import merge_metrics, resolve_csv_paths  # noqa: E402


def _write(d, name, df):
    df.to_csv(os.path.join(d, name), index=False)


def test_id_union_across_disjoint_files():
    """Same column, disjoint IDs across files -> union of rows, no dropping."""
    with tempfile.TemporaryDirectory() as d:
        _write(d, "1_c0_plddt.csv", pd.DataFrame({"ID": ["a", "b"], "mean_plddt": [95, 80]}))
        _write(d, "2_c1_plddt.csv", pd.DataFrame({"ID": ["c", "d"], "mean_plddt": [70, 60]}))
        df = merge_metrics([d]).set_index("ID")
    assert set(df.index) == {"a", "b", "c", "d"}
    assert df.loc["a", "mean_plddt"] == 95.0 and df.loc["d", "mean_plddt"] == 60.0


def test_first_wins_on_overlap():
    """A repeated column on overlapping IDs takes the first file's value per cell."""
    with tempfile.TemporaryDirectory() as d:
        _write(d, "1_first.csv", pd.DataFrame({"ID": ["a", "b"], "sequence": ["FA", "FB"]}))
        _write(d, "2_second.csv", pd.DataFrame({"ID": ["a", "b"], "sequence": ["SA", "SB"]}))
        df = merge_metrics([d]).set_index("ID")
    assert df.loc["a", "sequence"] == "FA" and df.loc["b", "sequence"] == "FB"


def test_first_wins_fills_gaps_not_overwrites():
    """First-wins keeps earlier values but still fills a later file's new rows/cells."""
    with tempfile.TemporaryDirectory() as d:
        # File 1 has a value for 'a' only; file 2 supplies 'b' (gap-fill) but must NOT
        # overwrite 'a'.
        _write(d, "1_x.csv", pd.DataFrame({"ID": ["a", "b"], "m": [1.0, np.nan]}))
        _write(d, "2_x.csv", pd.DataFrame({"ID": ["a", "b"], "m": [9.0, 2.0]}))
        df = merge_metrics([d]).set_index("ID")
    assert df.loc["a", "m"] == 1.0  # earlier file wins
    assert df.loc["b", "m"] == 2.0  # gap filled from later file


def test_id_alias_becomes_canonical_id():
    with tempfile.TemporaryDirectory() as d:
        _write(d, "ee.csv", pd.DataFrame({"id": ["a", "b"], "FPP_score": [0.9, 0.2]}))
        df = merge_metrics([d])
    assert "ID" in df.columns and "id" not in df.columns
    assert set(df["ID"]) == {"a", "b"}


def test_self_suffix():
    with tempfile.TemporaryDirectory() as d:
        _write(d, "g_max_sequence_identity.csv",
               pd.DataFrame({"ID": ["a"], "max_sequence_identity": [0.5]}))
        _write(d, "g_max_sequence_identity_self.csv",
               pd.DataFrame({"ID": ["a"], "max_sequence_identity": [0.9]}))
        df = merge_metrics([d]).set_index("ID")
    assert "max_sequence_identity" in df.columns
    assert "max_sequence_identity_self" in df.columns
    assert df.loc["a", "max_sequence_identity"] == 0.5
    assert df.loc["a", "max_sequence_identity_self"] == 0.9


def test_numeric_coercion_all_or_nothing():
    with tempfile.TemporaryDirectory() as d:
        _write(d, "x.csv", pd.DataFrame({
            "ID": ["a", "b"],
            "num": ["1.5", "2.5"],           # all numeric -> float
            "seq": ["MAA", "MBB"],           # strings -> stay object
            "flag": ["True", "False"],       # bool-ish strings -> NOT numeric, stay str
        }))
        df = merge_metrics([d])
    assert df["num"].dtype.kind == "f"
    assert df["seq"].dtype == object
    assert df["flag"].dtype == object


def test_missing_tokens_normalised():
    with tempfile.TemporaryDirectory() as d:
        _write(d, "x.csv", pd.DataFrame({"ID": ["a", "b", "c"], "m": ["1.0", "NA", "None"]}))
        df = merge_metrics([d]).set_index("ID")
    assert df.loc["a", "m"] == 1.0
    assert pd.isna(df.loc["b", "m"]) and pd.isna(df.loc["c", "m"])


def test_raw_feature_matrix_skipped():
    with tempfile.TemporaryDirectory() as d:
        wide = {"ID": ["a", "b"]}
        for i in range(300):
            wide[f"e{i}"] = [float(i), float(i)]
        _write(d, "emb.csv", pd.DataFrame(wide))
        _write(d, "plddt.csv", pd.DataFrame({"ID": ["a", "b"], "mean_plddt": [90, 80]}))
        df = merge_metrics([d])
    # The 300-col matrix is skipped; only the real metric survives.
    assert "mean_plddt" in df.columns and "e0" not in df.columns


def test_resolve_csv_paths_dedup_and_order():
    with tempfile.TemporaryDirectory() as d:
        _write(d, "a.csv", pd.DataFrame({"ID": ["x"], "m": [1]}))
        _write(d, "b.csv", pd.DataFrame({"ID": ["y"], "m": [2]}))
        one = os.path.join(d, "a.csv")
        paths = resolve_csv_paths([d, one])  # dir + explicit dup of a.csv
        assert len([p for p in paths if p.endswith("a.csv")]) == 1  # deduped
        assert all(p.endswith(".csv") for p in paths)


def test_no_usable_csv_raises():
    with tempfile.TemporaryDirectory() as d:
        try:
            merge_metrics([d])
        except ValueError as e:
            assert "no usable CSVs" in str(e)
        else:
            raise AssertionError("expected ValueError on empty dir")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
