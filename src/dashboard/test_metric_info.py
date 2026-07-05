from __future__ import annotations

"""Self-contained tests for dashboard/metric_info.py (per-metric metadata table).

Run from this directory (flat-module import, like build_dashboard does):
    cd src/dashboard && python test_metric_info.py
or under pytest:
    cd src/dashboard && python -m pytest test_metric_info.py -q

Pure metadata-consistency checks (no rendering): every metric's category is a
declared category; the metric-info and metric-category tables cover the same
metrics; and each info entry is well-formed (non-empty explanation + a columns
dict of string range-chips). "Other" must exist as the catch-all category.
"""

from metric_info import CATEGORY_ORDER, METRIC_CATEGORY, METRIC_INFO


def test_category_order_unique_and_has_other():
    assert len(CATEGORY_ORDER) == len(set(CATEGORY_ORDER)), "duplicate categories"
    assert "Other" in CATEGORY_ORDER, "catch-all 'Other' category missing"


def test_every_metric_category_value_is_declared():
    bad = {m: c for m, c in METRIC_CATEGORY.items() if c not in CATEGORY_ORDER}
    assert not bad, f"metrics assigned to undeclared categories: {bad}"


def test_info_and_category_cover_same_metrics():
    only_info = set(METRIC_INFO) - set(METRIC_CATEGORY)
    only_cat = set(METRIC_CATEGORY) - set(METRIC_INFO)
    assert not only_info, f"in METRIC_INFO but no category: {only_info}"
    assert not only_cat, f"has category but no METRIC_INFO entry: {only_cat}"


def test_info_entries_wellformed():
    for metric, entry in METRIC_INFO.items():
        assert isinstance(entry.get("explanation"), str) and entry["explanation"].strip(), metric
        cols = entry.get("columns")
        assert isinstance(cols, dict) and cols, f"{metric}: missing/empty columns"
        for cname, chip in cols.items():
            assert isinstance(cname, str) and cname, metric
            assert isinstance(chip, str), f"{metric}.{cname}: range chip not a string"


def test_no_duplicate_columns_within_metric():
    # dict keys are unique by construction; assert the count is sane (> 0) and
    # that column names don't carry stray whitespace.
    for metric, entry in METRIC_INFO.items():
        for cname in entry["columns"]:
            assert cname == cname.strip(), f"{metric}: '{cname}' has whitespace"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
