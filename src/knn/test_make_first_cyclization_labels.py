"""Unit tests for make_first_cyclization_labels.py (label-file construction).

Run: python test_make_first_cyclization_labels.py   (no pytest dependency required).

Uses only synthetic in-memory source tables. Exercises the per-enzyme collapse of
multi-product rows to ONE coarse label = most frequent first-cyclization class, with
ties broken toward the smallest class id (the determinism guarantee the k-NN metric
relies on).
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import make_first_cyclization_labels as mod  # noqa: E402


def _run(source_rows) -> pd.DataFrame:
    tmp = tempfile.mkdtemp(prefix="fc_labels_")
    src = os.path.join(tmp, "src.csv")
    out = os.path.join(tmp, "out.csv")
    pd.DataFrame(source_rows).to_csv(src, index=False)
    argv = sys.argv
    sys.argv = ["make_first_cyclization_labels.py", "--source", src, "--output", out]
    try:
        mod.main()
    finally:
        sys.argv = argv
    return pd.read_csv(out)


def test_most_frequent_class_wins():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 0},
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 0},
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 1},
        {"Enzyme_marts_ID": "marts_E2", "First_cyclization_product_id": 5},
    ]).set_index("reference_id")
    assert list(df.columns) == ["label"], df.columns
    assert int(df.loc["marts_E1", "label"]) == 0     # 0 appears twice, 1 once
    assert int(df.loc["marts_E2", "label"]) == 5
    print("ok most_frequent_class_wins")


def test_tie_breaks_to_smallest_class_id():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 7},
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 3},
    ]).set_index("reference_id")
    # 7 and 3 each appear once -> tie -> smallest id 3 wins (deterministic).
    assert int(df.loc["marts_E1", "label"]) == 3
    print("ok tie_breaks_to_smallest_class_id")


def test_nan_rows_dropped_and_int_cast():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 2.0},
        {"Enzyme_marts_ID": "marts_E2", "First_cyclization_product_id": None},
    ]).set_index("reference_id")
    assert int(df.loc["marts_E1", "label"]) == 2
    assert "marts_E2" not in df.index      # NaN label row dropped entirely
    print("ok nan_rows_dropped_and_int_cast")


def test_output_sorted_by_reference_id():
    df = _run([
        {"Enzyme_marts_ID": "marts_E3", "First_cyclization_product_id": 1},
        {"Enzyme_marts_ID": "marts_E1", "First_cyclization_product_id": 1},
        {"Enzyme_marts_ID": "marts_E2", "First_cyclization_product_id": 1},
    ])
    assert list(df["reference_id"]) == ["marts_E1", "marts_E2", "marts_E3"]
    print("ok output_sorted_by_reference_id")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
