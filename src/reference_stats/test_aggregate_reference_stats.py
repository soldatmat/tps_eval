from __future__ import annotations

"""Self-contained tests for aggregate_reference_stats.py.

Run from this directory (flat-module imports resolve like the runner does):
    cd src/reference_stats && python test_aggregate_reference_stats.py
or under pytest:
    cd src/reference_stats && python -m pytest test_aggregate_reference_stats.py -q

Synthetic in-memory Series / tiny temp CSVs whose summary statistics are known in
closed form, so the numeric stats (ddof=1 std, numpy-linear percentiles), the
numeric-vs-categorical classification, the per-class ``by_<labeling>`` stratification,
the CSV discovery (skip labeling + comparative files), the filename->metric mapping,
the label-file loader, and the NaN/Inf JSON-sanitiser are all checked exactly.
"""

import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import aggregate_reference_stats as agg  # noqa: E402


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b}"


def test_numeric_stats_closed_form():
    s = agg._numeric_stats(pd.Series([10.0, 20.0, 30.0, 40.0]))
    assert s["kind"] == "numeric"
    assert s["count"] == 4 and s["n_missing"] == 0
    _approx(s["mean"], 25.0)
    _approx(s["std"], math.sqrt(500.0 / 3.0))   # ddof=1 sample std
    _approx(s["min"], 10.0)
    _approx(s["median"], 25.0)
    _approx(s["max"], 40.0)
    _approx(s["p25"], 17.5)   # numpy linear interpolation
    _approx(s["p75"], 32.5)


def test_numeric_stats_missing_and_single():
    empty = agg._numeric_stats(pd.Series([np.nan, np.nan]))
    assert empty["count"] == 0 and empty["n_missing"] == 2
    assert empty["mean"] is None and empty["p25"] is None
    single = agg._numeric_stats(pd.Series([7.0]))
    assert single["count"] == 1
    _approx(single["std"], 0.0)   # ddof=1 undefined -> 0.0 by convention


def test_categorical_stats():
    s = agg._categorical_stats(pd.Series(["alpha", "alpha", "beta", np.nan]))
    assert s["kind"] == "categorical"
    assert s["count"] == 3 and s["n_missing"] == 1
    assert s["frequencies"] == {"alpha": 2, "beta": 1}
    assert s["n_unique"] == 2


def test_classify_bool_and_numeric_and_labels():
    assert agg._classify_and_stat(pd.Series([True, False, True]))["kind"] == "categorical"
    assert agg._classify_and_stat(pd.Series([1.0, 2.0, 3.0]))["kind"] == "numeric"
    # object column of numeric-looking strings -> numeric.
    assert agg._classify_and_stat(pd.Series(["1", "2", "3"]))["kind"] == "numeric"
    # object column of category labels -> categorical.
    assert agg._classify_and_stat(pd.Series(["single", "two"]))["kind"] == "categorical"


def test_normalize_id_column():
    df = pd.DataFrame({"id": ["a"], "m": [1]})
    assert "ID" in agg._normalize_id_column(df).columns
    # Already-canonical ID left untouched.
    df2 = pd.DataFrame({"ID": ["a"], "id": ["b"]})
    assert list(agg._normalize_id_column(df2).columns) == ["ID", "id"]


def test_aggregate_csv_drops_helpers_and_groups():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "structs_esmfold_radius_of_gyration.csv")
        pd.DataFrame({
            "ID": ["r1", "r2", "r3", "r4"],
            "radius_of_gyration": [10.0, 20.0, 30.0, 40.0],
            "sequence": ["A", "B", "C", "D"],   # in _DROP_COLUMNS -> excluded
        }).to_csv(p, index=False)
        result = agg.aggregate_csv(p, labelings={
            "arch": {"r1": "single", "r2": "single", "r3": "two", "r4": "two"}})
    cols = result["columns"]
    assert "radius_of_gyration" in cols
    assert "sequence" not in cols and "ID" not in cols
    rg = cols["radius_of_gyration"]
    _approx(rg["by_arch"]["single"]["mean"], 15.0)
    _approx(rg["by_arch"]["two"]["mean"], 35.0)
    assert result["n_rows"] == 4


def test_discover_csvs_skips_labels_and_comparative():
    with tempfile.TemporaryDirectory() as d:
        pd.DataFrame({"ID": ["a"], "radius_of_gyration": [1.0]}).to_csv(
            os.path.join(d, "structs_radius_of_gyration.csv"), index=False)
        pd.DataFrame({"reference_id": ["a"], "label": ["x"]}).to_csv(
            os.path.join(d, "arch_labels.csv"), index=False)
        pd.DataFrame({"ID": ["a"], "max_identity": [0.9]}).to_csv(
            os.path.join(d, "g_swissprot_search.csv"), index=False)
        found = agg.discover_csvs(d)
    assert "radius_of_gyration" in found
    assert not any("swissprot" in m for m in found)
    assert not any("labels" in os.path.basename(p) for p in found.values())


def test_metric_from_filename():
    assert agg._metric_from_filename("g_motifs.csv") == ("motif_search", False)
    assert agg._metric_from_filename("structs_esmfold_plddt.csv") == ("plddt", False)
    m, comp = agg._metric_from_filename("gen_swissprot_search.csv")
    assert comp is True and m == "swissprot_search"
    # Unknown metric: a known input prefix is stripped so it still bands.
    m2, comp2 = agg._metric_from_filename("structs_my_new_metric.csv")
    assert comp2 is False and m2 == "my_new_metric"


def test_load_label_file():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "arch_labels.csv")
        pd.DataFrame({"reference_id": ["r1", "r2", "r3", "r2"],
                      "label": ["single", "", "two", "double"]}).to_csv(p, index=False)
        mapping = agg.load_label_file(p)
    assert mapping["r1"] == "single"
    assert "r2" in mapping and mapping["r2"] == "double"   # blank dropped, last dup wins
    assert mapping["r3"] == "two"


def test_load_label_file_too_few_columns_raises():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "bad.csv")
        pd.DataFrame({"only": ["r1"]}).to_csv(p, index=False)
        try:
            agg.load_label_file(p)
        except SystemExit:
            pass
        else:
            raise AssertionError("expected SystemExit for <2-column label file")


def test_labeling_name_from_path():
    assert agg.labeling_name_from_path("/x/first_cyclization_labels.csv") == "first_cyclization"
    assert agg.labeling_name_from_path("/x/domain.csv") == "domain"


def test_json_safe():
    safe = agg._json_safe({"a": float("nan"), "b": [1.0, float("inf")], "c": 2.0})
    assert safe["a"] is None
    assert safe["b"] == [1.0, None]
    assert safe["c"] == 2.0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
