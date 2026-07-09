"""Unit tests for make_substrate_labels.py (substrate-class label-file construction).

Run: python test_make_substrate_labels.py   (no pytest dependency required).

Synthetic in-memory MARTS `Type` tables only. Exercises the Type->substrate mapping,
lower/strip normalization, dropping of un-mapped Type values, and the per-enzyme
collapse to ONE label = most frequent class with ties broken toward the SMALLEST
carbon count (the determinism guarantee).
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd


import tps_eval.knn.make_substrate_labels as mod  # noqa: E402
from tps_eval.knn.make_substrate_labels import TYPE_TO_SUBSTRATE, _carbons  # noqa: E402


def _run(source_rows) -> pd.DataFrame:
    tmp = tempfile.mkdtemp(prefix="sub_labels_")
    src = os.path.join(tmp, "src.csv")
    out = os.path.join(tmp, "out.csv")
    pd.DataFrame(source_rows).to_csv(src, index=False)
    argv = sys.argv
    sys.argv = ["make_substrate_labels.py", "--source", src, "--output", out]
    try:
        mod.main()
    finally:
        sys.argv = argv
    return pd.read_csv(out)


def test_type_mapping_and_carbons():
    assert TYPE_TO_SUBSTRATE["mono"] == "GPP"
    assert TYPE_TO_SUBSTRATE["sesq"] == "FPP"
    assert TYPE_TO_SUBSTRATE["di"] == "GGPP"
    assert TYPE_TO_SUBSTRATE["sqs"] == "EDSQ"     # squalene synthase folds to EDSQ
    assert _carbons("GPP") < _carbons("FPP") < _carbons("GGPP")
    assert _carbons("not_a_class") == 1000        # unknown sorts last
    print("ok type_mapping_and_carbons")


def test_most_frequent_class_wins():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "Type": "mono"},
        {"Enzyme_marts_ID": "marts_E1", "Type": "mono"},
        {"Enzyme_marts_ID": "marts_E1", "Type": "sesq"},
    ]).set_index("reference_id")
    assert df.loc["marts_E1", "label"] == "GPP"   # mono x2 beats sesq x1
    print("ok most_frequent_class_wins")


def test_tie_breaks_to_smallest_carbon():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "Type": "di"},    # GGPP, C20
        {"Enzyme_marts_ID": "marts_E1", "Type": "tri"},   # EDSQ, C30
    ]).set_index("reference_id")
    # 1:1 tie -> smallest carbon count wins: GGPP (20) over EDSQ (30).
    assert df.loc["marts_E1", "label"] == "GGPP"
    print("ok tie_breaks_to_smallest_carbon")


def test_case_and_whitespace_normalized():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "Type": "  SESQ "},
    ]).set_index("reference_id")
    assert df.loc["marts_E1", "label"] == "FPP"
    print("ok case_and_whitespace_normalized")


def test_unmapped_type_dropped():
    df = _run([
        {"Enzyme_marts_ID": "marts_E1", "Type": "mono"},
        {"Enzyme_marts_ID": "marts_E2", "Type": "bogus_type"},
    ]).set_index("reference_id")
    assert "marts_E1" in df.index and df.loc["marts_E1", "label"] == "GPP"
    assert "marts_E2" not in df.index      # only-unmapped enzyme drops out entirely
    print("ok unmapped_type_dropped")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
