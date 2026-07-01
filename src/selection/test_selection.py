"""Self-contained tests for the selection primitives (gate / score / band_filter /
diversity_dedup) and the merge util.

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_selection.py
or under pytest:
    cd src/selection && python -m pytest test_selection.py -q

The diversity_dedup test is skipped when the ``mmseqs`` binary is not on PATH.
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from band_filter import apply_band_filter
from gate import apply_gate
from merge import merge_metrics
from score import apply_score


def _df():
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "class": ["c0", "c0", "c1", "c1"],
        "mean_plddt": [95.0, 80.0, 92.0, np.nan],
        "proteinmpnn_nll": [1.5, 2.5, 1.8, 1.6],
        "af3_iptm": [0.95, 0.90, 0.93, 0.94],
        "mg_canonical_motif_coordination": [True, False, True, True],
        "diphosphate_to_nearest_ion": [2.0, 2.5, 5.0, 2.8],
        "sequence_identity": [0.5, 0.6, 0.71, 0.8],
        "domain_architecture": ["two", "two", "single", "single"],
        "catalytic_pocket_volume": [800, 200, 1000, 1500],
    })


def test_gate_and_and_missing():
    out, rep = apply_gate(_df(), [
        {"col": "mg_canonical_motif_coordination", "eq": True},
        {"col": "diphosphate_to_nearest_ion", "le": 3.0},
    ])
    assert set(out["ID"]) == {"a", "d"}, set(out["ID"])
    assert rep["n_pass"] == 2
    print("ok gate AND + boolean + le")


def test_gate_novelty_threshold():
    out, _ = apply_gate(_df(), [{"col": "sequence_identity", "le": 0.72}])
    assert set(out["ID"]) == {"a", "b", "c"}, set(out["ID"])
    print("ok gate novelty threshold")


def test_gate_missing_fails():
    out, _ = apply_gate(_df(), [{"col": "mean_plddt", "ge": 90}])
    assert set(out["ID"]) == {"a", "c"}, set(out["ID"])  # d NaN fails; b 80 fails
    print("ok gate missing-value fails")


def test_gate_when_conditional():
    df = _df()
    # c10-only-style: enforce catalytic_pocket_volume>=900 ONLY for class c1; c0 rows unaffected.
    out, _ = apply_gate(df, [{"when": {"col": "class", "eq": "c1"},
                              "col": "catalytic_pocket_volume", "ge": 1200}])
    # c0 rows (a,b) pass regardless; c1 rows: c(1000)<1200 fails, d(1500) passes.
    assert set(out["ID"]) == {"a", "b", "d"}, set(out["ID"])
    print("ok gate when-conditional (class-specific)")


def test_gate_any_of():
    out, _ = apply_gate(_df(), [{"any_of": [
        {"col": "mean_plddt", "ge": 94}, {"col": "af3_iptm", "ge": 0.935}]}])
    assert set(out["ID"]) == {"a", "d"}, set(out["ID"])
    print("ok gate nested any_of")


def test_score_within_group_and_direction():
    df = _df().dropna(subset=["mean_plddt"]).reset_index(drop=True)
    out, _ = apply_score(df, [
        {"col": "mean_plddt", "weight": 1, "direction": "higher"},
        {"col": "proteinmpnn_nll", "weight": 1, "direction": "lower"},
    ], zscore_within="class")
    c0 = out[out["class"] == "c0"].set_index("ID")
    assert c0.loc["a", "score"] > c0.loc["b", "score"]
    assert c0.loc["a", "score_rank"] == 1
    print("ok score z-sum within group + direction flip")


def test_band_filter_per_architecture():
    out, _ = apply_band_filter(_df(), {
        "catalytic_pocket_volume": {"by": "domain_architecture",
                                    "single": {"min": 617, "max": 1377},
                                    "two": {"min": 326, "max": 1016}},
    })
    assert set(out["ID"]) == {"a", "c"}, set(out["ID"])
    print("ok band_filter per-architecture")


def test_band_filter_leaf_onesided():
    out, _ = apply_band_filter(_df(), {"af3_iptm": {"min": 0.93}})
    assert set(out["ID"]) == {"a", "c", "d"}, set(out["ID"])
    print("ok band_filter one-sided leaf")


def test_merge_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        pd.DataFrame({"ID": ["a", "b"], "mean_plddt": [95, 80]}).to_csv(
            os.path.join(d, "g_plddt.csv"), index=False)
        pd.DataFrame({"id": ["a", "b"], "FPP_score": [0.9, 0.2], "sequence": ["MAA", "MBB"]}
                     ).to_csv(os.path.join(d, "g_enzyme_explorer_sequence_only.csv"), index=False)
        df = merge_metrics([d])
    assert set(df.columns) >= {"ID", "mean_plddt", "FPP_score", "sequence"}
    assert len(df) == 2
    assert df.set_index("ID").loc["a", "mean_plddt"] == 95.0
    print("ok merge round-trip (id->ID, numeric coerce, sequence kept)")


def test_export_bands_overall_and_by():
    from export_bands import export_bands
    ref = {"reference_set": "marts_db", "structure_source": "esmfold", "metrics": {
        "pocket_descriptors": {"columns": {
            "catalytic_pocket_volume": {"kind": "numeric", "p25": 300.0, "p75": 1000.0,
                "by_arch": {"single": {"kind": "numeric", "p25": 617.0, "p75": 1377.0},
                            "two": {"kind": "numeric", "p25": 326.0, "p75": 1016.0}}}}}}}
    overall = export_bands(ref, {"catalytic_pocket_volume": {}}, {"lo": "p25", "hi": "p75"})
    assert overall["metrics"]["catalytic_pocket_volume"] == {"min": 300.0, "max": 1000.0}
    byarch = export_bands(ref, {"catalytic_pocket_volume": {}}, {"lo": "p25", "hi": "p75"}, by="arch")
    leaf = byarch["metrics"]["catalytic_pocket_volume"]
    assert leaf["by"] == "arch" and leaf["single"] == {"min": 617.0, "max": 1377.0}
    print("ok export_bands overall + per-architecture (band_filter format)")


def test_diversity_dedup_per_group():
    if shutil.which("mmseqs") is None:
        print("SKIP diversity_dedup (mmseqs not on PATH)")
        return
    from diversity_dedup import apply_diversity_dedup
    df = pd.DataFrame({
        "ID": ["a1", "a2", "b1"],
        "class": ["c0", "c0", "c1"],
        "sequence": ["MKKAAAILLVVGG", "MKKAAAILLVVGG", "WQWQWQWQWQWQ"],
        "score": [1.0, 0.5, 0.9],
    })
    out, _ = apply_diversity_dedup(df, quality_col="score", id_threshold=0.9,
                                   group_col="class")
    assert "a1" in set(out["ID"]) and "a2" not in set(out["ID"]), set(out["ID"])
    assert "b1" in set(out["ID"])
    print("ok diversity_dedup best-rep-per-cluster per group")


if __name__ == "__main__":
    test_gate_and_and_missing()
    test_gate_novelty_threshold()
    test_gate_missing_fails()
    test_gate_when_conditional()
    test_gate_any_of()
    test_score_within_group_and_direction()
    test_band_filter_per_architecture()
    test_band_filter_leaf_onesided()
    test_merge_roundtrip()
    test_export_bands_overall_and_by()
    test_diversity_dedup_per_group()
    print("ALL TESTS PASSED")
