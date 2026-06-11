"""Unit tests for the label-agnostic k-NN coarse-label transfer logic.

Run: python test_knn_label_transfer.py   (no pytest dependency required).
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from knn_label_transfer import (  # noqa: E402
    ABSTAIN_LABEL,
    _strip_chain_suffix,
    calibrate,
    score_to_similarity,
    transfer_labels,
    vote_space,
)


def _write_topk(path, rows):
    pd.DataFrame(rows, columns=["query_id", "rank", "neighbour_id", "score"]).to_csv(
        path, index=False
    )


def test_score_to_similarity():
    assert abs(score_to_similarity("sequence", 40.0) - 0.40) < 1e-9
    assert abs(score_to_similarity("structural", 0.5) - 0.5) < 1e-9
    assert abs(score_to_similarity("embedding", 0.0) - 1.0) < 1e-9
    assert abs(score_to_similarity("embedding", 1.0) - 0.5) < 1e-9
    print("ok score_to_similarity")


def test_chain_suffix():
    valid = {"marts_E00123", "weird_id"}
    assert _strip_chain_suffix("marts_E00123_A", valid) == "marts_E00123"
    assert _strip_chain_suffix("marts_E00123", valid) == "marts_E00123"
    # full id is itself valid (contains underscore) -> not stripped
    assert _strip_chain_suffix("weird_id", valid) == "weird_id"
    # no valid set: strip short trailing token
    assert _strip_chain_suffix("foo_A", None) == "foo"
    assert _strip_chain_suffix("foo_bar", None) == "foo_bar"  # token too long
    print("ok chain_suffix")


def test_vote_unanimous_and_abstain():
    label_map = {"r1": "A", "r2": "A", "r3": "B"}
    classes = ["A", "B"]
    nbrs = [(1, "r1", 90.0), (2, "r2", 80.0), (3, "r3", 50.0)]
    v = vote_space(nbrs, "sequence", label_map, classes, tau=0.40)
    assert v.predicted == "A", v.predicted
    assert v.n_voters == 3
    assert 0 < v.confidence <= 1
    # all below tau -> abstain
    nbrs_low = [(1, "r1", 10.0), (2, "r3", 5.0)]
    v2 = vote_space(nbrs_low, "sequence", label_map, classes, tau=0.40)
    assert v2.predicted is None
    assert v2.n_voters == 0
    print("ok vote unanimous + abstain")


def test_structural_chain_join():
    label_map = {"marts_E1": "A", "marts_E2": "B"}
    classes = ["A", "B"]
    valid = set(label_map)
    nbrs = [(1, "marts_E1_A", 0.9), (2, "marts_E2_B", 0.8)]
    v = vote_space(nbrs, "structural", label_map, classes, tau=0.5, valid_ids=valid)
    assert v.predicted == "A"  # higher TM-score -> heavier vote
    assert v.n_voters == 2
    print("ok structural chain join")


def test_calibrate_and_predict_end_to_end():
    tmp = tempfile.mkdtemp(prefix="knn_test_")
    # 6 references, 2 perfectly-separated classes; neighbours within class are close.
    label_file = os.path.join(tmp, "labels.csv")
    pd.DataFrame(
        {"reference_id": [f"r{i}" for i in range(6)],
         "label": ["A", "A", "A", "B", "B", "B"]}
    ).to_csv(label_file, index=False)

    # SELF top-k: each ref's neighbours are its same-class refs (high sim) then
    # cross-class (low sim). sequence space (identity %).
    seq_self = os.path.join(tmp, "seq_self_topk.csv")
    rows = []
    same = {"A": ["r0", "r1", "r2"], "B": ["r3", "r4", "r5"]}
    lbl = {"r0": "A", "r1": "A", "r2": "A", "r3": "B", "r4": "B", "r5": "B"}
    for q in [f"r{i}" for i in range(6)]:
        cls = lbl[q]
        rank = 1
        for n in same[cls]:
            if n == q:
                continue
            rows.append({"query_id": q, "rank": rank, "neighbour_id": n, "score": 85.0})
            rank += 1
        # one far cross-class neighbour
        other = "r3" if cls == "A" else "r0"
        rows.append({"query_id": q, "rank": rank, "neighbour_id": other, "score": 20.0})
    _write_topk(seq_self, rows)

    cal = calibrate({"sequence": seq_self}, label_file, labeling="test")
    assert cal["spaces"]["sequence"]["accuracy"] == 1.0, cal["spaces"]["sequence"]
    assert cal["ensemble"]["accuracy"] == 1.0

    # PREDICT: a design close to class A, and a novel design far from everything.
    seq_design = os.path.join(tmp, "seq_design_topk.csv")
    _write_topk(seq_design, [
        {"query_id": "d_A", "rank": 1, "neighbour_id": "r0", "score": 80.0},
        {"query_id": "d_A", "rank": 2, "neighbour_id": "r1", "score": 78.0},
        {"query_id": "d_novel", "rank": 1, "neighbour_id": "r0", "score": 12.0},
        {"query_id": "d_novel", "rank": 2, "neighbour_id": "r3", "score": 10.0},
    ])
    df = transfer_labels({"sequence": seq_design}, label_file, cal)
    df = df.set_index("ID")
    assert df.loc["d_A", "predicted_label"] == "A", df.loc["d_A"].to_dict()
    assert df.loc["d_novel", "predicted_label"] == ABSTAIN_LABEL, df.loc["d_novel"].to_dict()
    assert df.loc["d_novel", "confidence"] == 0.0
    print("ok calibrate + predict end-to-end (incl. abstain)")


if __name__ == "__main__":
    test_score_to_similarity()
    test_chain_suffix()
    test_vote_unanimous_and_abstain()
    test_structural_chain_join()
    test_calibrate_and_predict_end_to_end()
    print("ALL TESTS PASSED")
