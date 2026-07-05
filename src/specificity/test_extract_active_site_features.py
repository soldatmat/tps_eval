from __future__ import annotations

"""Unit tests for the active-site feature-extraction logic (extract_active_site_features).

Run: python test_extract_active_site_features.py   (no pytest dependency required).

Builds tiny synthetic .pdb structures with a real DDXXD-family motif and explicit
side-chain coordinating oxygens, so the metal-point / shell-selection / property-
profile path runs end-to-end without real MARTS-DB data. Also tests the alignment-
free composition math (closed-form fractions), the no-motif NaN contract (OSC
outlier), ID keying, CSV filename, and empty/malformed -> NaN row (not a crash).
"""

import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from extract_active_site_features import (  # noqa: E402
    AA_ORDER,
    FEATURE_COLUMNS,
    PROFILE_COLUMNS,
    _default_save_path,
    _nan_features,
    active_site_features,
    active_site_features_dir,
)

AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}


def _atom_line(serial: int, atom_name: str, resname: str, resseq: int, xyz, element: str) -> str:
    x, y, z = xyz
    region = " " + f"{atom_name:^4}" + " "  # cols 12-17; cols 13-16 hold the name
    return (
        f"ATOM  {serial:5d}{region}{resname} A{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element}"
    )


def _write_pdb(path: str, seq: str, ca_coords, asp_oxygens=None) -> None:
    """One CA per residue; for aspartate residues in `asp_oxygens` also write OD1/OD2
    (the metal-coordinating side-chain oxygens the tool anchors the metal point on).

    asp_oxygens: {residue_index: ((od1_xyz), (od2_xyz))}.
    """
    asp_oxygens = asp_oxygens or {}
    lines = []
    serial = 1
    for i, (aa, ca) in enumerate(zip(seq, ca_coords)):
        lines.append(_atom_line(serial, "CA", AA3[aa], i + 1, ca, "C"))
        serial += 1
        if i in asp_oxygens:
            od1, od2 = asp_oxygens[i]
            lines.append(_atom_line(serial, "OD1", AA3[aa], i + 1, od1, "O"))
            serial += 1
            lines.append(_atom_line(serial, "OD2", AA3[aa], i + 1, od2, "O"))
            serial += 1
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def test_nan_features():
    d = _nan_features()
    assert d["n_coordinating_oxygens"] == 0
    # Every OTHER feature column is NaN.
    for c in FEATURE_COLUMNS:
        if c == "n_coordinating_oxygens":
            continue
        assert math.isnan(d[c]), c
    print("ok _nan_features")


def test_default_save_path():
    assert _default_save_path("/x/gen_structs") == "/x/gen_structs_active_site_features.csv"
    print("ok _default_save_path")


def test_features_with_ddxxd_motif():
    """DDXXD present + coordinating oxygens -> metal_point_found, closed-form
    composition fractions over the shell, correct oxygen count + inner aromatics."""
    with tempfile.TemporaryDirectory() as d:
        # DDAAD matches [DE][DE]..[DE] at index 0-4; coordinating offsets (0,1,4)
        # -> residues 0,1,4 (all D). Trailing FYWK gives aromatics + a basic residue.
        seq = "DDAADFYWK"  # D D A A D F Y W K
        # All CAs within 12 A of the (near-origin) metal point -> whole seq in shell.
        ca = [(float(i), 0.0, 0.0) for i in range(len(seq))]
        asp_ox = {
            0: ((0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
            1: ((1.0, 1.0, 0.0), (1.0, -1.0, 0.0)),
            4: ((4.0, 1.0, 0.0), (4.0, -1.0, 0.0)),
        }
        p = os.path.join(d, "des1.pdb")
        _write_pdb(p, seq, ca, asp_ox)

        feat = active_site_features(p, radius=12.0)
        assert feat["metal_point_found"] is True, feat
        assert feat["n_residues"] == 9, feat
        assert feat["n_shell_residues"] == 9, feat            # all within 12 A
        assert feat["n_coordinating_oxygens"] == 6, feat      # 3 Asp * OD1/OD2

        # Closed-form composition over the 9 shell residues (D3 A2 F1 Y1 W1 K1).
        assert abs(feat["frac_acidic"] - 3 / 9) < 1e-9, feat
        assert abs(feat["frac_aromatic"] - 3 / 9) < 1e-9, feat   # F,Y,W
        assert abs(feat["frac_aliphatic"] - 2 / 9) < 1e-9, feat  # A,A
        assert abs(feat["frac_basic"] - 1 / 9) < 1e-9, feat      # K
        assert abs(feat["frac_polar"]) < 1e-9, feat
        assert abs(feat["frac_glycine"]) < 1e-9, feat
        assert abs(feat["frac_proline"]) < 1e-9, feat
        # The mutually-exclusive property partition sums to 1 over the shell.
        assert abs(sum(feat[c] for c in PROFILE_COLUMNS) - 1.0) < 1e-9, feat
        # Per-AA fractions.
        assert abs(feat["frac_aa_D"] - 3 / 9) < 1e-9, feat
        assert abs(feat["frac_aa_A"] - 2 / 9) < 1e-9, feat
        # The 20 per-AA fractions sum to 1 (no non-standard residues here).
        assert abs(sum(feat[f"frac_aa_{a}"] for a in AA_ORDER) - 1.0) < 1e-9, feat

        # Inner-shell aromatics: F,Y,W CA are at x=5,6,7 -> within 8 A of the
        # ~origin metal point; K (basic) excluded.
        assert feat["n_aromatic_within_8A"] == 3, feat
        # Geometry descriptors are finite and positive.
        assert feat["carboxylate_convergence_radius"] > 0, feat
        assert feat["shell_radius_of_gyration"] > 0, feat
        assert feat["mean_dist_to_metal_point"] > 0, feat
    print("ok features with DDXXD motif (composition + geometry)")


def test_shell_radius_restricts():
    """A tight radius keeps only residues near the metal point in the shell."""
    with tempfile.TemporaryDirectory() as d:
        seq = "DDAADAAAAAAAAAAAAAAA"  # DDXXD at 0-4, then a long tail marching away
        ca = [(float(i) * 3.0, 0.0, 0.0) for i in range(len(seq))]
        asp_ox = {
            0: ((0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
            1: ((3.0, 1.0, 0.0), (3.0, -1.0, 0.0)),
            4: ((12.0, 1.0, 0.0), (12.0, -1.0, 0.0)),
        }
        p = os.path.join(d, "des2.pdb")
        _write_pdb(p, seq, ca, asp_ox)
        wide = active_site_features(p, radius=12.0)
        narrow = active_site_features(p, radius=5.0)
        assert narrow["n_shell_residues"] < wide["n_shell_residues"], (narrow, wide)
        assert narrow["metal_point_found"] is True
    print("ok shell radius restricts membership")


def test_no_motif_nan_row():
    """No DDXXD -> no metal point -> metal_point_found False, features NaN (OSC outlier)."""
    with tempfile.TemporaryDirectory() as d:
        seq = "AAAAAAAAAA"
        ca = [(float(i) * 3.8, 0.0, 0.0) for i in range(len(seq))]
        p = os.path.join(d, "osc.pdb")
        _write_pdb(p, seq, ca)
        feat = active_site_features(p, radius=12.0)
        assert feat["metal_point_found"] is False, feat
        assert feat["n_shell_residues"] == 0, feat
        for c in FEATURE_COLUMNS:
            if c == "n_coordinating_oxygens":
                assert feat[c] == 0, feat
            else:
                assert math.isnan(feat[c]), (c, feat[c])
    print("ok no-motif -> NaN row")


def test_dir_driver_keys_and_csv():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "gen_structs")
        os.makedirs(structs)
        # One good (motif) structure, one OSC (no motif).
        seq_good = "DDAADFYWK"
        ca_good = [(float(i), 0.0, 0.0) for i in range(len(seq_good))]
        _write_pdb(os.path.join(structs, "good.pdb"), seq_good, ca_good, {
            0: ((0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
            1: ((1.0, 1.0, 0.0), (1.0, -1.0, 0.0)),
            4: ((4.0, 1.0, 0.0), (4.0, -1.0, 0.0)),
        })
        _write_pdb(os.path.join(structs, "osc.pdb"), "AAAAAAAA",
                   [(float(i), 0.0, 0.0) for i in range(8)])
        # A malformed structure -> parse failure -> NaN row, not a crash.
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("this is not a valid pdb file\n")

        df = active_site_features_dir(structs)
        assert list(df["id"]) == sorted(df["id"]), "rows sorted by id"
        assert set(df["id"]) == {"good", "osc", "broken"}, df["id"].tolist()
        by_id = df.set_index("id")
        assert bool(by_id.loc["good", "metal_point_found"]) is True
        assert bool(by_id.loc["osc", "metal_point_found"]) is False
        assert bool(by_id.loc["broken", "metal_point_found"]) is False
        assert pd.isna(by_id.loc["osc", "frac_aromatic"])
        # CSV written to the default sibling path.
        assert os.path.exists(_default_save_path(structs))
    print("ok dir driver: id keying + CSV + malformed -> NaN")


def main():
    test_nan_features()
    test_default_save_path()
    test_features_with_ddxxd_motif()
    test_shell_radius_restricts()
    test_no_motif_nan_row()
    test_dir_driver_keys_and_csv()
    print("\nAll 6 tests passed.")


if __name__ == "__main__":
    main()
