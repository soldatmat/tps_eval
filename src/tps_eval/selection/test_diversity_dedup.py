from __future__ import annotations

"""Self-contained tests for diversity_dedup.py — WITHOUT the real mmseqs binary.

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_diversity_dedup.py
or under pytest:
    cd src/selection && python -m pytest test_diversity_dedup.py -q

The mmseqs shell-out is the only impure part; we never run it. Instead we either
(a) monkeypatch ``subprocess.run`` to capture the command line and drop a synthetic
``*_cluster.tsv`` (so the command construction + TSV parsing in ``_mmseqs_clusters``
is exercised), or (b) monkeypatch ``_mmseqs_clusters`` directly to test the pure
grouping / best-rep / top-N / threshold logic. No mmseqs is invoked.
"""

import os
import sys
import types

import numpy as np
import pandas as pd


import tps_eval.selection.diversity_dedup as dd  # noqa: E402
from tps_eval.selection.diversity_dedup import apply_diversity_dedup  # noqa: E402


def _fake_run_factory(clusters):
    """Return a subprocess.run stand-in that writes ``clusters`` (list of (rep, member))
    to ``<out_prefix>_cluster.tsv`` and records the command line."""
    recorded = {}

    def fake_run(cmd, capture_output, text):
        recorded["cmd"] = cmd
        # cmd = [mmseqs, easy-cluster, <fasta>, <out_prefix>, <tmp>, ...]
        out_prefix = cmd[3]
        with open(out_prefix + "_cluster.tsv", "w") as fh:
            for rep, member in clusters:
                fh.write(f"{rep}\t{member}\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    return fake_run, recorded


def test_command_construction_and_tsv_parse(monkeypatch=None):
    # a1,a2 cluster together (rep a1); b1 alone.
    fake_run, recorded = _fake_run_factory([("a1", "a1"), ("a1", "a2"), ("b1", "b1")])
    orig = dd.subprocess.run
    dd.subprocess.run = fake_run
    try:
        df = pd.DataFrame({"ID": ["a1", "a2", "b1"],
                           "sequence": ["MK", "MK", "WQ"], "score": [1.0, 0.5, 0.9]})
        out, _ = apply_diversity_dedup(df, quality_col="score", id_threshold=0.9,
                                       coverage=0.75)
    finally:
        dd.subprocess.run = orig
    cmd = recorded["cmd"]
    assert cmd[0] == "mmseqs" and cmd[1] == "easy-cluster"
    assert cmd[cmd.index("--min-seq-id") + 1] == "0.9"
    assert cmd[cmd.index("-c") + 1] == "0.75"
    assert "--threads" in cmd
    # a1 (best of its cluster) + b1 survive; a2 (lower score) is dropped.
    assert set(out["ID"]) == {"a1", "b1"}


def test_best_rep_per_cluster(monkeypatch=None):
    def fake_clusters(id_to_seq, min_seq_id, coverage=0.8, threads=4):
        # all three in one cluster represented by x1
        return {k: "x1" for k in id_to_seq}
    orig = dd._mmseqs_clusters
    dd._mmseqs_clusters = fake_clusters
    try:
        df = pd.DataFrame({"ID": ["x1", "x2", "x3"],
                           "sequence": ["A", "A", "A"], "q": [0.2, 0.9, 0.5]})
        out, rep = apply_diversity_dedup(df, quality_col="q", id_threshold=0.9)
    finally:
        dd._mmseqs_clusters = orig
    # single cluster -> only the highest-quality member (x2) survives.
    assert set(out["ID"]) == {"x2"}
    assert rep["n_in"] == 3 and rep["n_out"] == 1


def test_per_group_threshold_and_top_n(monkeypatch=None):
    def fake_clusters(id_to_seq, min_seq_id, coverage=0.8, threads=4):
        return {k: k for k in id_to_seq}  # every seq its own cluster
    orig = dd._mmseqs_clusters
    dd._mmseqs_clusters = fake_clusters
    try:
        df = pd.DataFrame({
            "ID": ["a1", "a2", "a3", "b1"],
            "cls": ["c0", "c0", "c0", "c1"],
            "sequence": ["AA", "BB", "CC", "DD"],
            "q": [0.3, 0.9, 0.6, 0.7],
        })
        out, rep = apply_diversity_dedup(
            df, quality_col="q", group_col="cls",
            id_threshold_per_group={"c0": 0.5, "c1": 0.8},
            n_out_per_group=2)
    finally:
        dd._mmseqs_clusters = orig
    # No clustering collapses (all own clusters), so top-2 by q per group:
    # c0 -> a2(0.9), a3(0.6); c1 -> b1.
    assert set(out[out["cls"] == "c0"]["ID"]) == {"a2", "a3"}
    assert set(out[out["cls"] == "c1"]["ID"]) == {"b1"}
    groups = {g["group"]: g for g in rep["groups"]}
    assert groups["c0"]["min_seq_id"] == 0.5 and groups["c1"]["min_seq_id"] == 0.8


def test_singleton_group_no_mmseqs():
    # A single-sequence group must NOT shell out to mmseqs (guarded by len<=1).
    def boom(*a, **k):
        raise AssertionError("mmseqs must not be called for a singleton group")
    orig = dd._mmseqs_clusters
    dd._mmseqs_clusters = boom
    try:
        df = pd.DataFrame({"ID": ["only"], "sequence": ["MK"], "q": [1.0]})
        out, _ = apply_diversity_dedup(df, quality_col="q", id_threshold=0.9)
    finally:
        dd._mmseqs_clusters = orig
    assert set(out["ID"]) == {"only"}


def test_missing_threshold_raises():
    df = pd.DataFrame({"ID": ["a", "b"], "sequence": ["A", "B"], "q": [1.0, 0.5]})
    try:
        apply_diversity_dedup(df, quality_col="q")  # no threshold at all
    except ValueError as e:
        assert "id_threshold" in str(e)
    else:
        raise AssertionError("expected ValueError when no threshold given")


def test_missing_sequence_column_raises():
    df = pd.DataFrame({"ID": ["a"], "q": [1.0]})
    try:
        apply_diversity_dedup(df, quality_col="q", id_threshold=0.9)
    except ValueError as e:
        assert "sequence" in str(e)
    else:
        raise AssertionError("expected ValueError for missing sequence column")


def test_missing_quality_column_raises():
    df = pd.DataFrame({"ID": ["a"], "sequence": ["MK"]})
    try:
        apply_diversity_dedup(df, quality_col="q", id_threshold=0.9)
    except ValueError as e:
        assert "quality_col" in str(e)
    else:
        raise AssertionError("expected ValueError for missing quality column")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
