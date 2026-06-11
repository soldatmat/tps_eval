from __future__ import annotations

"""Self-contained tests for domain_structural_identity.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_domain_structural_identity.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_domain_structural_identity.py -q

These exercise the pure-logic surface that does NOT need EnzymeExplorer or foldseek
installed: the module_id <-> (id, type) parsing and the per-hit -> per-design
reduction (best TM-score, matched reference + type, per-type bests, NaN rows for
zero-domain designs). The end-to-end detect+align path is validated separately on a
compute node (see scripts/<cluster>/jobs/domain_structural_identity.sh).
"""

import os
import tempfile

import numpy as np
import pandas as pd

from domain_structural_identity import (
    COLUMNS,
    DOMAIN_TYPES,
    _parse_module_id,
    _reduce_alignments,
)


def test_parse_module_id_basic():
    assert _parse_module_id("mydesign_alpha_0") == ("mydesign", "alpha")
    assert _parse_module_id("marts_E00002_beta_0") == ("marts_E00002", "beta")
    # IDs with underscores keep their underscores; type is the 2nd-to-last field.
    assert _parse_module_id("run_0001_sample_3_gamma_2") == ("run_0001_sample_3", "gamma")


def test_parse_module_id_all_types():
    for t in DOMAIN_TYPES:
        assert _parse_module_id(f"design7_{t}_0") == ("design7", t)


def test_parse_module_id_unrecognised():
    # No trailing _<type>_<index> -> whole stem is the id, empty type.
    assert _parse_module_id("just_a_name") == ("just_a_name", "")
    # Looks like a number suffix but type token is not a known domain type.
    assert _parse_module_id("foo_widget_0") == ("foo_widget_0", "")


def _write_hits(path: str, rows):
    cols = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart",
            "qend", "tstart", "tend", "evalue", "bits", "alntmscore", "qtmscore",
            "ttmscore", "lddt"]
    df = pd.DataFrame([{**{c: 0 for c in cols}, **r} for r in rows], columns=cols)
    df.to_csv(path, index=False)


def test_reduce_best_and_per_type():
    tmp = tempfile.mkdtemp()
    raw = os.path.join(tmp, "hits.csv")
    # design "d1" has an alpha domain matching two alpha refs and a beta ref.
    _write_hits(raw, [
        {"query": "d1_alpha_0", "target": "marts_E1_alpha_0", "alntmscore": 0.55, "lddt": 0.60},
        {"query": "d1_alpha_0", "target": "marts_E2_alpha_0", "alntmscore": 0.81, "lddt": 0.77},
        {"query": "d1_beta_0", "target": "marts_E3_beta_0", "alntmscore": 0.40, "lddt": 0.50},
    ])
    counts = {"d1": 2, "d2_nodomain": 0}
    df = _reduce_alignments(raw, counts)

    assert list(df.columns) == COLUMNS
    assert set(df["ID"]) == {"d1", "d2_nodomain"}

    r1 = df[df["ID"] == "d1"].iloc[0]
    # Global best is the 0.81 alpha hit to marts_E2_alpha_0.
    assert abs(r1["domain_structural_tmscore_to_known"] - 0.81) < 1e-9
    assert r1["domain_structural_tmscore_to_known_hit"] == "marts_E2_alpha_0"
    assert r1["domain_structural_tmscore_to_known_type"] == "alpha"
    assert abs(r1["domain_structural_lddt_to_known"] - 0.77) < 1e-9
    assert r1["n_detected_domains"] == 2
    # Per-type bests.
    assert abs(r1["domain_structural_tmscore_to_known_alpha"] - 0.81) < 1e-9
    assert abs(r1["domain_structural_tmscore_to_known_beta"] - 0.40) < 1e-9
    assert np.isnan(r1["domain_structural_tmscore_to_known_gamma"])

    # Zero-domain design -> NaN metrics, n_detected_domains == 0.
    r2 = df[df["ID"] == "d2_nodomain"].iloc[0]
    assert r2["n_detected_domains"] == 0
    assert np.isnan(r2["domain_structural_tmscore_to_known"])
    assert pd.isna(r2["domain_structural_tmscore_to_known_hit"])


def test_reduce_design_with_no_hits_is_nan():
    # A design that was detected (count>0) but produced NO foldseek hits -> NaN.
    tmp = tempfile.mkdtemp()
    raw = os.path.join(tmp, "hits.csv")
    _write_hits(raw, [
        {"query": "other_alpha_0", "target": "marts_E1_alpha_0", "alntmscore": 0.7, "lddt": 0.6},
    ])
    counts = {"detected_but_unmatched": 1, "other": 1}
    df = _reduce_alignments(raw, counts)
    r = df[df["ID"] == "detected_but_unmatched"].iloc[0]
    assert r["n_detected_domains"] == 1
    assert np.isnan(r["domain_structural_tmscore_to_known"])


if __name__ == "__main__":
    test_parse_module_id_basic()
    test_parse_module_id_all_types()
    test_parse_module_id_unrecognised()
    test_reduce_best_and_per_type()
    test_reduce_design_with_no_hits_is_nan()
    print("All domain_structural_identity tests passed.")
