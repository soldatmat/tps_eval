"""Unit tests for extract_ee_domain_features.py (EE domain-comparison feature block).

Run: python test_extract_ee_domain_features.py   (no EnzymeExplorer / Foldseek needed).

Tests the pure matrix-assembly logic with synthetic mocks: ID sanitization (the
underscore-in-ID bug guard), the reference-module column union across fold classifiers,
and the per-protein ``1 - best-TM-score`` fill including the alpha1/alpha2 detected-domain
split and the same-type-only matching rule. The EE detection/comparison subprocess calls
are NOT exercised.
"""
from __future__ import annotations

import os
import sys
from types import SimpleNamespace


from tps_eval.enzyme_explorer.extract_ee_domain_features import (
    _sanitize_id,
    build_feature_matrix,
    build_module_columns,
)


def test_sanitize_id():
    assert _sanitize_id("marts_E00000") == "martsE00000"
    assert _sanitize_id("a_b-c.d") == "abcd"
    assert _sanitize_id("clean123") == "clean123"
    print("ok sanitize_id")


def test_build_module_columns_union_and_order():
    clf1 = SimpleNamespace(domain_type_2_order_of_domain_modules={
        "alpha1": [("mA", 0), ("mB", 1)],
        "beta": [("mC", 0)],
    })
    clf2 = SimpleNamespace(domain_type_2_order_of_domain_modules={
        "alpha1": [("mB", 0), ("mZ", 1)],   # mB shared; mZ new
        "beta": [("mC", 0)],
    })
    ordered, module_to_type = build_module_columns([clf1, clf2])
    # Union of alpha1 modules (sorted) precedes beta modules (DOMAIN_TYPE_ORDER).
    assert ordered == ["mA", "mB", "mZ", "mC"], ordered
    assert module_to_type["mA"] == "alpha1"
    assert module_to_type["mC"] == "beta"
    print("ok build_module_columns_union_and_order")


def _det(domain, module_id):
    return SimpleNamespace(domain=domain, module_id=module_id)


def test_build_feature_matrix_one_minus_tm_and_alpha_split():
    module_columns = ["mA1", "mA2", "mB"]
    module_to_type = {"mA1": "alpha1", "mA2": "alpha2", "mB": "beta"}
    san_to_orig = {"martsE1": "marts_E1"}
    detected_domains = {
        "martsE1": [
            _det("alpha", "martsE1_alpha_0"),   # first alpha -> alpha1
            _det("alpha", "martsE1_alpha_1"),   # second alpha -> alpha2
            _det("beta", "martsE1_beta_0"),
        ]
    }
    comparison_results = {
        "martsE1": {
            "martsE1_alpha_0": [("mA1", 0.8), ("mB", 0.99)],  # mB is cross-type -> ignored
            "martsE1_alpha_1": [("mA2", 0.6)],
            "martsE1_beta_0": [("mB", 0.7)],
        }
    }
    df = build_feature_matrix(
        ["marts_E1"], san_to_orig, detected_domains, comparison_results,
        module_columns, module_to_type,
    ).set_index("id")
    # value = 1 - best same-type TM.
    assert abs(df.loc["marts_E1", "mA1"] - (1 - 0.8)) < 1e-6
    assert abs(df.loc["marts_E1", "mA2"] - (1 - 0.6)) < 1e-6
    assert abs(df.loc["marts_E1", "mB"] - (1 - 0.7)) < 1e-6  # NOT 1-0.99 (cross-type ignored)
    print("ok build_feature_matrix_one_minus_tm_and_alpha_split")


def test_build_feature_matrix_no_comparison_gives_ones():
    module_columns = ["mA1"]
    module_to_type = {"mA1": "alpha1"}
    df = build_feature_matrix(
        ["marts_E9"], {"martsE9": "marts_E9"}, {}, {}, module_columns, module_to_type,
    ).set_index("id")
    # No comparison -> feat stays 0 -> emitted 1 - 0 = 1.0.
    assert abs(df.loc["marts_E9", "mA1"] - 1.0) < 1e-6
    print("ok build_feature_matrix_no_comparison_gives_ones")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
