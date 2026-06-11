"""Unit tests for the SDR-divergence logic.

Run: python test_sdr_divergence.py   (no pytest dependency required).

Builds tiny synthetic .pdb structures so the structural-superposition + position-
mapping path is exercised end-to-end without the real MARTS-DB / design data.
"""
from __future__ import annotations

import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sdr_divergence import (  # noqa: E402
    Panel,
    ResidueInfo,
    _rank1_neighbours,
    _to_similarity,
    load_panel,
    sdr_divergence_one,
)

AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}


def _write_pdb(path: str, seq: str, coords) -> None:
    """One CA atom per residue at the given coords (enough for superposition + the
    nearest-CA mapping; the metal-point/structure-derived panel path is tested
    separately with the explicit panel here so we don't need side chains)."""
    lines = []
    for i, (aa, xyz) in enumerate(zip(seq, coords), start=1):
        x, y, z = xyz
        lines.append(
            f"ATOM  {i:5d}  CA  {AA3[aa]} A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def test_to_similarity():
    assert abs(_to_similarity("sequence", 60.0) - 0.6) < 1e-9
    assert abs(_to_similarity("structural", 0.72) - 0.72) < 1e-9
    assert math.isnan(_to_similarity("structural", float("nan")))
    print("ok _to_similarity")


def test_panel_match_and_indices():
    panel = Panel({"5eat": [(10, "T"), (12, "Y"), (99, "A")]})
    assert panel._match_reference("5eat") == "5eat"
    assert panel._match_reference("5eat_A") == "5eat"  # chain-suffixed neighbour
    assert panel._match_reference("totally_other") is None
    # build a structure whose author resnums are 1..15
    with tempfile.TemporaryDirectory() as d:
        seq = "ACDEFGHIKLMNPQR"
        coords = [(float(i), 0.0, 0.0) for i in range(len(seq))]
        p = os.path.join(d, "5eat.pdb")
        _write_pdb(p, seq, coords)
        info = ResidueInfo(p)
        idx = panel.indices_for("5eat", info)
        # resnums 10 and 12 are present (1-based authoring), 99 is not.
        assert idx == [9, 11], idx
    print("ok panel match + indices")


def test_divergence_flag_and_identity():
    """Design identical in backbone to the neighbour but with substitutions at two of
    three panel positions -> sdr_identity = 1/3, and (high global sim) -> flagged."""
    with tempfile.TemporaryDirectory() as d:
        nbr_seq = "MTTYVKLANDE"   # neighbour
        # design: same length/backbone, differs at panel positions 2 (T->S) and 4 (Y->F)
        des_seq = "MSTFVKLANDE"
        coords = [(float(i) * 3.8, 0.0, 0.0) for i in range(len(nbr_seq))]
        np_path = os.path.join(d, "neighbour.pdb")
        dp_path = os.path.join(d, "design.pdb")
        _write_pdb(np_path, nbr_seq, coords)
        _write_pdb(dp_path, des_seq, coords)
        nbr = ResidueInfo(np_path)
        des = ResidueInfo(dp_path)
        # explicit panel over resnums 2,3,4 (1-based author numbering)
        panel = Panel({"neighbour": [(2, "T"), (3, "T"), (4, "Y")]})
        row = sdr_divergence_one(
            des, "neighbour", nbr, similarity=0.78, space="structural",
            panel=panel, tau_high=0.6, tau_low=0.7,
        )
        assert row["n_sdr_positions"] == 3, row
        assert abs(row["sdr_identity"] - (1.0 / 3.0)) < 1e-9, row
        assert row["n_sdr_mismatches"] == 2, row
        assert row["specificity_divergence"] is True, row
        # divergent_positions records neighbour-residue + resnum + design-residue
        assert "T2S" in row["divergent_positions"], row["divergent_positions"]
        assert "Y4F" in row["divergent_positions"], row["divergent_positions"]
        print("ok divergence flag + identity:", row["divergent_positions"])

        # Same design but global sim below tau_high -> NOT flagged (wrong regime).
        row2 = sdr_divergence_one(
            des, "neighbour", nbr, similarity=0.50, space="structural",
            panel=panel, tau_high=0.6, tau_low=0.7,
        )
        # sdr_identity still computed, but flag must be False.
        assert row2["specificity_divergence"] is False, row2
        print("ok below-tau_high not flagged")

        # Identical design (perfect SDR match) -> high identity, NOT flagged.
        row3 = sdr_divergence_one(
            nbr, "neighbour", nbr, similarity=0.95, space="structural",
            panel=panel, tau_high=0.6, tau_low=0.7,
        )
        assert abs(row3["sdr_identity"] - 1.0) < 1e-9, row3
        assert row3["specificity_divergence"] is False, row3
        print("ok identical design not flagged")


def test_structure_derived_panel_needs_motifs():
    """A backbone-only synthetic structure has no DDXXD/NSE motifs in sequence, so the
    structure-derived panel returns None -> NaN row (n_sdr_positions == 0)."""
    with tempfile.TemporaryDirectory() as d:
        seq = "AAAAAAAAAA"
        coords = [(float(i) * 3.8, 0.0, 0.0) for i in range(len(seq))]
        p = os.path.join(d, "x.pdb")
        _write_pdb(p, seq, coords)
        info = ResidueInfo(p)
        row = sdr_divergence_one(
            info, "x", info, similarity=0.9, space="structural", panel=None,
        )
        assert row["n_sdr_positions"] == 0, row
        assert math.isnan(row["sdr_identity"]), row
        assert row["specificity_divergence"] is False, row
    print("ok structure-derived panel NaN without motifs")


def test_rank1_neighbours_chain_strip_and_preference():
    with tempfile.TemporaryDirectory() as d:
        seq_csv = os.path.join(d, "seq_topk.csv")
        str_csv = os.path.join(d, "str_topk.csv")
        pd.DataFrame(
            [["q1", 1, "marts_E1", 55.0], ["q1", 2, "marts_E2", 40.0]],
            columns=["query_id", "rank", "neighbour_id", "score"],
        ).to_csv(seq_csv, index=False)
        pd.DataFrame(
            [["q1", 1, "marts_E9_A", 0.82], ["q2", 1, "marts_E3_B", 0.71]],
            columns=["query_id", "rank", "neighbour_id", "score"],
        ).to_csv(str_csv, index=False)
        valid = {"marts_E1", "marts_E2", "marts_E3", "marts_E9"}
        nn = _rank1_neighbours(str_csv, seq_csv, valid)
        # q1: structural preferred over sequence, chain suffix stripped
        assert nn["q1"] == ("marts_E9", 0.82, "structural"), nn["q1"]
        # q2: only structural, chain suffix stripped
        assert nn["q2"] == ("marts_E3", 0.71, "structural"), nn["q2"]
    print("ok rank1 neighbours: chain-strip + structural preference")


def test_load_panel_skips_comments():
    """The committed starter panel (with '#' provenance lines) loads cleanly."""
    here = os.path.dirname(os.path.abspath(__file__))
    starter = os.path.join(here, "sdr_panel_teas_hps.csv")
    if os.path.exists(starter):
        panel = load_panel(starter)
        assert "5eat" in panel.references(), panel.references()
        assert len(panel.by_ref["5eat"]) == 9, panel.by_ref["5eat"]
        print("ok starter panel loads (9 positions, comments skipped)")
    else:
        print("skip starter-panel test (file absent)")


if __name__ == "__main__":
    test_to_similarity()
    test_panel_match_and_indices()
    test_divergence_flag_and_identity()
    test_structure_derived_panel_needs_motifs()
    test_rank1_neighbours_chain_strip_and_preference()
    test_load_panel_skips_comments()
    print("\nAll SDR-divergence tests passed.")
