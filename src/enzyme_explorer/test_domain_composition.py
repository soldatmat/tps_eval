"""Unit tests for domain_composition.py (TPS domain composition table builder).

Run: python test_domain_composition.py   (no pytest / EnzymeExplorer needed).

Tests only the pure, EE-free logic: regions_to_rows (per-type counts, ordered
architecture string, zero-domain fill, unexpected-type handling), the JSON sidecar
parser, and the default save-path derivation. Detection itself (detect_domains_json)
imports EnzymeExplorer lazily and is NOT exercised here.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from domain_composition import (  # noqa: E402
    COLUMNS,
    _default_save_path,
    load_detections_json,
    regions_to_rows,
)


def test_regions_to_rows_counts_and_architecture():
    seq_to_regions = {
        "s1": [
            {"domain": "beta", "module_id": "s1_beta_0"},
            {"domain": "alpha", "module_id": "s1_alpha_1"},
            {"domain": "alpha", "module_id": "s1_alpha_0"},
        ],
    }
    df = regions_to_rows(seq_to_regions, ["s1", "s2"]).set_index("ID")
    assert list(df.reset_index().columns) == COLUMNS
    # s1: 3 domains, 2 alpha, 1 beta; architecture sorted by (type, index).
    assert int(df.loc["s1", "n_domains"]) == 3
    assert int(df.loc["s1", "n_alpha"]) == 2
    assert int(df.loc["s1", "n_beta"]) == 1
    assert df.loc["s1", "domain_architecture"] == "alpha-alpha-beta"
    # s2: absent from detections -> zero-domain row.
    assert int(df.loc["s2", "n_domains"]) == 0
    assert df.loc["s2", "domain_architecture"] == ""
    print("ok regions_to_rows_counts_and_architecture")


def test_unexpected_domain_type_counted_without_column():
    seq_to_regions = {"s1": [{"domain": "weird", "module_id": "s1_weird_0"},
                             {"domain": "alpha", "module_id": "s1_alpha_0"}]}
    df = regions_to_rows(seq_to_regions, ["s1"]).set_index("ID")
    # Unexpected type still counts toward n_domains + architecture but has no n_* column.
    assert int(df.loc["s1", "n_domains"]) == 2
    assert "weird" in df.loc["s1", "domain_architecture"]
    assert "n_weird" not in df.columns
    print("ok unexpected_domain_type_counted_without_column")


def test_load_detections_json_roundtrip():
    tmp = tempfile.mkdtemp(prefix="domcomp_")
    p = os.path.join(tmp, "det.json")
    payload = {"s1": [{"module_id": "s1_alpha_0", "domain": "alpha", "tmscore": 0.9}]}
    with open(p, "w") as fh:
        json.dump(payload, fh)
    parsed = load_detections_json(p)
    assert parsed["s1"][0]["domain"] == "alpha"
    print("ok load_detections_json_roundtrip")


def test_default_save_path():
    assert _default_save_path("/a/b/structs").endswith("structs_domain_composition.csv")
    # Trailing separator handled.
    assert _default_save_path("/a/b/structs/").endswith("structs_domain_composition.csv")
    print("ok default_save_path")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
