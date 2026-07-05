from __future__ import annotations

"""Self-contained tests for select_designs.py (the composite JSON-spec driver).

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_select_designs.py
or under pytest:
    cd src/selection && python -m pytest test_select_designs.py -q

Synthetic in-memory DataFrames + temp output dir. Checks JSON-spec composition
(gate -> score -> take_top_n survivors + reports), group_from_id regex synthesis,
the per-group top-N cap by score, sequence injection from a fasta_map, the written
survivors.csv / survivors.fasta / provenance manifest.md, the diversity_dedup op
wiring (mmseqs mocked away), and the unknown-op guard.
"""

import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import select_designs as sd  # noqa: E402
from select_designs import run_selection, select_and_write, write_manifest  # noqa: E402


def _df():
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d", "e"],
        "class": ["c0", "c0", "c1", "c1", "c1"],
        "mean_plddt": [95.0, 80.0, 92.0, 88.0, 91.0],
        "nll": [1.5, 2.5, 1.8, 1.6, 1.7],
        "sequence": ["MAAA", "MBBB", "MCCC", "MDDD", "MEEE"],
    })


def test_composition_gate_then_score_then_cap():
    spec = {
        "group_by": "class",
        "n_out_per_group": 1,
        "ops": [
            {"op": "gate", "conditions": [{"col": "mean_plddt", "ge": 90}]},
            {"op": "score", "terms": [{"col": "mean_plddt", "weight": 1, "direction": "higher"}]},
        ],
    }
    surv, reports, gb = run_selection(_df(), spec)
    assert gb == "class"
    # gate keeps a(95),c(92),e(91); top-1/class by plddt -> c0:a, c1:c.
    assert set(surv["ID"]) == {"a", "c"}
    assert [r["op"] for r in reports] == ["input", "gate", "score", "take_top_n"]
    gate_rep = reports[1]
    assert gate_rep["op"] == "gate" and gate_rep["n_pass"] == 3


def test_group_from_id_regex():
    df = pd.DataFrame({
        "ID": ["gen_c0_1", "gen_c0_2", "gen_c1_1"],
        "mean_plddt": [95.0, 80.0, 92.0],
        "sequence": ["A", "B", "C"],
    })
    spec = {"group_by": "class", "group_from_id": r"_(c\d+)_", "n_out_per_group": 1,
            "ops": [{"op": "score", "terms": [{"col": "mean_plddt", "direction": "higher"}]}]}
    surv, _, gb = run_selection(df, spec)
    assert gb == "class"
    assert set(surv["class"]) == {"c0", "c1"}
    # top-1 per synthesised class: c0 -> gen_c0_1 (95), c1 -> gen_c1_1.
    assert set(surv["ID"]) == {"gen_c0_1", "gen_c1_1"}


def test_sequence_injected_from_fasta_map():
    df = pd.DataFrame({"ID": ["a", "b"], "mean_plddt": [95.0, 80.0]})  # no sequence col
    spec = {"n_out_per_group": 2, "ops": []}
    surv, _, _ = run_selection(df, spec, fasta_map={"a": "MAAA", "b": "MBBB"})
    assert "sequence" in surv.columns
    assert set(surv["sequence"]) == {"MAAA", "MBBB"}


def test_take_top_n_by_score_no_group():
    df = pd.DataFrame({"ID": ["a", "b", "c"], "mean_plddt": [70.0, 90.0, 80.0],
                       "sequence": ["A", "B", "C"]})
    spec = {"n_out_per_group": 2,
            "ops": [{"op": "score", "terms": [{"col": "mean_plddt", "direction": "higher"}]}]}
    surv, _, _ = run_selection(df, spec)
    # global top-2 by plddt: b(90), c(80).
    assert set(surv["ID"]) == {"b", "c"}


def test_diversity_dedup_op_wiring():
    # mmseqs is mocked away; verify select routes the op with group_col defaulting to
    # group_by and passes quality_col / id_threshold through.
    captured = {}

    def fake_dedup(df, **kwargs):
        captured.update(kwargs)
        return df, {"op": "diversity_dedup", "n_in": len(df), "n_out": len(df), "groups": []}

    orig = sd.apply_diversity_dedup
    sd.apply_diversity_dedup = fake_dedup
    try:
        spec = {"group_by": "class",
                "ops": [{"op": "diversity_dedup", "quality_col": "mean_plddt",
                         "id_threshold": 0.7}]}
        run_selection(_df(), spec)
    finally:
        sd.apply_diversity_dedup = orig
    assert captured["quality_col"] == "mean_plddt"
    assert captured["id_threshold"] == 0.7
    assert captured["group_col"] == "class"  # defaulted to group_by


def test_unknown_op_raises():
    try:
        run_selection(_df(), {"ops": [{"op": "frobnicate"}]})
    except ValueError as e:
        assert "unknown selection op" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown op")


def test_select_and_write_outputs_and_manifest():
    spec = {
        "group_by": "class",
        "n_out_per_group": 1,
        "ops": [
            {"op": "gate", "conditions": [{"col": "mean_plddt", "ge": 90}]},
            {"op": "score", "terms": [
                {"col": "mean_plddt", "weight": 1, "direction": "higher"},
                {"col": "nll", "weight": 1, "direction": "lower"}]},
        ],
    }
    with tempfile.TemporaryDirectory() as d:
        prefix = os.path.join(d, "sel")
        surv = select_and_write(_df(), spec, prefix, title="Phase-test")
        assert os.path.isfile(prefix + "_survivors.csv")
        assert os.path.isfile(prefix + "_survivors.fasta")
        assert os.path.isfile(prefix + "_manifest.md")
        # FASTA holds exactly the survivors' sequences, keyed by ID.
        fasta = open(prefix + "_survivors.fasta").read()
        for rid in surv["ID"]:
            assert f">{rid}\n" in fasta
        manifest = open(prefix + "_manifest.md").read()
        # Provenance records the ops + the gate condition + the score formula.
        assert "gate" in manifest and "score" in manifest
        assert "mean_plddt ge 90" in manifest
        assert "z(mean_plddt)" in manifest


def test_write_manifest_op_rows():
    reports = [
        {"op": "input", "n_in": 10, "group_counts": {"all": 10}},
        {"op": "gate", "n_in": 10, "n_pass": 6,
         "conditions": [{"condition": "x ge 1", "passed": 6}], "group_counts": {"all": 6}},
    ]
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "m.md")
        write_manifest(reports, {"ops": []}, out)
        text = open(out).read()
    assert "provenance manifest" in text
    assert "x ge 1 → 6" in text


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
