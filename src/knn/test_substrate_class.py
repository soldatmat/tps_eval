"""Unit tests for the substrate-class combiner.

Run: python test_substrate_class.py   (no pytest dependency required).
"""
from __future__ import annotations

import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from knn_label_transfer import calibrate  # noqa: E402
from substrate_class import (  # noqa: E402
    _within_one_size_class,
    combine_substrate_class,
    load_ee_substrate,
    pocket_volume_to_band,
)


def test_pocket_volume_band():
    assert pocket_volume_to_band(200.0) == "GPP"      # small -> C10
    assert pocket_volume_to_band(700.0) == "FPP"      # mid -> C15
    assert pocket_volume_to_band(1200.0) == "GGPP"    # larger -> C20
    assert pocket_volume_to_band(5000.0) == "GFPP"    # largest band
    assert pocket_volume_to_band(float("nan")) is None
    assert pocket_volume_to_band(None) is None
    print("ok pocket_volume_band")


def test_within_one_size_class():
    assert _within_one_size_class("GPP", "FPP") is True    # C10 vs C15 adjacent
    assert _within_one_size_class("GPP", "GGPP") is False   # C10 vs C20 (2 apart)
    assert _within_one_size_class("FPP", "FPP") is True
    assert _within_one_size_class("FPP", "IDS") is None      # IDS has no size rank
    print("ok within_one_size_class")


def test_ee_column_folding():
    # EE seq-only console-script schema: short *_score codes + TPS gate.
    tmp = tempfile.mkdtemp(prefix="subclass_ee_")
    ee_csv = os.path.join(tmp, "ee.csv")
    pd.DataFrame([
        # design d1: GGPP-ish (CPP folds into GGPP and dominates)
        {"ID": "d1", "FPP_score": 0.1, "GPP_score": 0.05, "GGPP_score": 0.2,
         "CPP_score": 0.6, "2xFPP_score": 0.01, "EDSQ_score": 0.02,
         "TPS_score": 0.99, "isTPS": True},
        # design d2: FPP dominates
        {"ID": "d2", "FPP_score": 0.7, "GPP_score": 0.1, "GGPP_score": 0.05,
         "CPP_score": 0.02, "2xFPP_score": 0.01, "EDSQ_score": 0.02,
         "TPS_score": 0.99, "isTPS": True},
    ]).to_csv(ee_csv, index=False)
    ee = load_ee_substrate(ee_csv)
    assert ee["d1"][0] == "GGPP", ee["d1"]   # CPP (0.6) folded into GGPP, beats raw GGPP 0.2
    assert ee["d2"][0] == "FPP", ee["d2"]
    # TPS gate must NOT be selected as a substrate
    assert ee["d1"][0] not in ("TPS", "isTPS")
    print("ok ee_column_folding")


def test_combiner_end_to_end():
    tmp = tempfile.mkdtemp(prefix="subclass_")
    # substrate label file: 6 refs, 2 classes (GPP small, GGPP large).
    label_file = os.path.join(tmp, "labels.csv")
    pd.DataFrame(
        {"reference_id": [f"r{i}" for i in range(6)],
         "label": ["GPP", "GPP", "GPP", "GGPP", "GGPP", "GGPP"]}
    ).to_csv(label_file, index=False)
    lbl = {"r0": "GPP", "r1": "GPP", "r2": "GPP", "r3": "GGPP", "r4": "GGPP", "r5": "GGPP"}
    same = {"GPP": ["r0", "r1", "r2"], "GGPP": ["r3", "r4", "r5"]}

    # SELF top-k for calibration (sequence space, identity %).
    seq_self = os.path.join(tmp, "seq_self_topk.csv")
    rows = []
    for q in [f"r{i}" for i in range(6)]:
        cls = lbl[q]
        rank = 1
        for n in same[cls]:
            if n == q:
                continue
            rows.append({"query_id": q, "rank": rank, "neighbour_id": n, "score": 85.0})
            rank += 1
        other = "r3" if cls == "GPP" else "r0"
        rows.append({"query_id": q, "rank": rank, "neighbour_id": other, "score": 20.0})
    pd.DataFrame(rows).to_csv(seq_self, index=False)
    cal = calibrate({"sequence": seq_self}, label_file, labeling="substrate_class")

    # PREDICT design top-k: d_small close to GPP class.
    seq_design = os.path.join(tmp, "seq_design_topk.csv")
    pd.DataFrame([
        {"query_id": "d_small", "rank": 1, "neighbour_id": "r0", "score": 82.0},
        {"query_id": "d_small", "rank": 2, "neighbour_id": "r1", "score": 79.0},
    ]).to_csv(seq_design, index=False)

    # pocket: small volume -> GPP band (agrees with k-NN GPP).
    pocket_csv = os.path.join(tmp, "pocket.csv")
    pd.DataFrame([{"ID": "d_small", "catalytic_pocket_volume": 250.0}]).to_csv(
        pocket_csv, index=False)

    # EE: GPP argmax (agrees).
    ee_csv = os.path.join(tmp, "ee.csv")
    pd.DataFrame([{"ID": "d_small", "GPP_score": 0.8, "FPP_score": 0.1,
                   "GGPP_score": 0.05, "TPS_score": 0.99}]).to_csv(ee_csv, index=False)

    df = combine_substrate_class(
        {"sequence": seq_design}, label_file, cal,
        pocket_csv=pocket_csv, ee_csv=ee_csv,
    ).set_index("ID")
    r = df.loc["d_small"]
    assert r["predicted_substrate"] == "GPP", r.to_dict()
    assert r["predicted_substrate_source"] == "knn"
    assert r["knn_substrate"] == "GPP"
    assert r["pocket_volume_band"] == "GPP"
    assert r["substrate_agreement"] is True or r["substrate_agreement"] == True  # noqa: E712
    assert r["ee_substrate"] == "GPP"
    assert r["ee_agreement"] is True or r["ee_agreement"] == True  # noqa: E712
    assert int(r["n_signals_agree"]) == 2
    print("ok combiner end-to-end (k-NN + pocket band + EE agreement)")


def test_ee_fallback_when_knn_abstains():
    tmp = tempfile.mkdtemp(prefix="subclass_fb_")
    label_file = os.path.join(tmp, "labels.csv")
    pd.DataFrame({"reference_id": ["r0", "r1", "r2", "r3", "r4", "r5"],
                  "label": ["GPP", "GPP", "GPP", "FPP", "FPP", "FPP"]}).to_csv(
        label_file, index=False)
    lbl = {"r0": "GPP", "r1": "GPP", "r2": "GPP", "r3": "FPP", "r4": "FPP", "r5": "FPP"}
    same = {"GPP": ["r0", "r1", "r2"], "FPP": ["r3", "r4", "r5"]}
    seq_self = os.path.join(tmp, "seq_self_topk.csv")
    rows = []
    for q in [f"r{i}" for i in range(6)]:
        cls = lbl[q]
        rank = 1
        for n in same[cls]:
            if n == q:
                continue
            rows.append({"query_id": q, "rank": rank, "neighbour_id": n, "score": 90.0})
            rank += 1
    pd.DataFrame(rows).to_csv(seq_self, index=False)
    cal = calibrate({"sequence": seq_self}, label_file, labeling="substrate_class")

    # design is novel -> all neighbours below tau -> k-NN abstains.
    seq_design = os.path.join(tmp, "seq_design_topk.csv")
    pd.DataFrame([
        {"query_id": "d_novel", "rank": 1, "neighbour_id": "r0", "score": 8.0},
        {"query_id": "d_novel", "rank": 2, "neighbour_id": "r3", "score": 6.0},
    ]).to_csv(seq_design, index=False)
    ee_csv = os.path.join(tmp, "ee.csv")
    pd.DataFrame([{"ID": "d_novel", "FPP_score": 0.7, "GPP_score": 0.1,
                   "TPS_score": 0.9}]).to_csv(ee_csv, index=False)

    df = combine_substrate_class({"sequence": seq_design}, label_file, cal,
                                 ee_csv=ee_csv).set_index("ID")
    r = df.loc["d_novel"]
    assert r["knn_substrate"] == "unknown", r.to_dict()
    assert r["predicted_substrate"] == "FPP"            # EE fallback
    assert r["predicted_substrate_source"] == "enzyme_explorer"
    assert r["ee_agreement"] == "NA"                     # no k-NN call to compare against
    print("ok EE fallback when k-NN abstains")


if __name__ == "__main__":
    test_pocket_volume_band()
    test_within_one_size_class()
    test_ee_column_folding()
    test_combiner_end_to_end()
    test_ee_fallback_when_knn_abstains()
    print("ALL TESTS PASSED")
