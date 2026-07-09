from __future__ import annotations

"""Self-contained tests for diphosphate_sensor.py (numpy + biopython only — no conda
env / EnzymeExplorer / PyMOL). Run from this directory so the flat-module imports
resolve:
    cd src/structure_metrics && python test_diphosphate_sensor.py
or:
    cd src/structure_metrics && python -m pytest test_diphosphate_sensor.py -q

Builds a synthetic protein whose residue names spell one DDXXD + one NSE/DTE motif whose
coordinating side-chain oxygens are placed at 12 sphere points summing to zero, so the
carboxylate-cage metal point sits at the ORIGIN. Basic (Arg/Lys) and Tyr residues are
then planted at known coordinates relative to that origin to exercise: the near+points-
toward counting of diphosphate-sensor residues, the direction check (side chain pointing
away is rejected), the distance cutoff, and both the sequence-adjacency and spatial RY-
pair criteria. Also covers the metal-point-absent -> all-zero/False contract, ID keying,
the sibling-CSV filename, and NaN/0-on-broken-structure.
"""

import os
import tempfile

import numpy as np

from tps_eval.structure_metrics.diphosphate_sensor import (
    COLUMNS,
    DEFAULT_CUTOFF,
    diphosphate_sensor_dir,
    diphosphate_sensor_one,
)

# ResidueInfo is the parsed-structure view the tool consumes (from the specificity module).
import sys
from pathlib import Path
from tps_eval.specificity.sdr_divergence import ResidueInfo  # noqa: E402

_THREE = {"A": "ALA", "D": "ASP", "N": "ASN", "L": "LEU", "S": "SER", "E": "GLU",
          "R": "ARG", "K": "LYS", "Y": "TYR"}

# Two-motif sequence (DDXXD at 5, NSE/DTE at 15) with R/K/Y planted in the leading pad.
_SEQ = "RKYAA" + "DDAAD" + "AAAAA" + "NDLASAAAE" + "AA"
_COORD_IDX = [5, 6, 9, 15, 16, 19, 23]
_OX_ATOMS = {"ASP": ("OD1", "OD2"), "ASN": ("OD1",), "SER": ("OG",), "GLU": ("OE1", "OE2")}

_R = 3.0
_AXIS6 = [(_R, 0, 0), (-_R, 0, 0), (0, _R, 0), (0, -_R, 0), (0, 0, _R), (0, 0, -_R)]
_SPHERE12 = _AXIS6 + _AXIS6


def _atom_line(serial, name, resname, chain, resseq, xyz, *, element=None):
    x, y, z = xyz
    if element is None:
        element = "O" if name.startswith("O") else ("N" if name.startswith("N") else name[0])
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{'ATOM':<6}{serial:>5} {atom_field}{' ':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}{'':10}{element:>2}\n"
    )


def _cage_oxygens_by_residue(seq=_SEQ):
    """{residue_index: [(atom_name, xyz)]} for the 12 coordinating oxygens (metal point
    at the origin)."""
    slots = []
    for i in _COORD_IDX:
        for name in _OX_ATOMS[_THREE[seq[i]]]:
            slots.append((i, name))
    assert len(slots) == len(_SPHERE12), (len(slots), len(_SPHERE12))
    out = {}
    for (i, name), xyz in zip(slots, _SPHERE12):
        out.setdefault(i, []).append((name, xyz))
    return out


def _write(path, seq, extra_atoms):
    """Write the cage + a CA per residue on a radius-20 circle (far from the origin so
    they don't interfere) + any ``extra_atoms`` ({idx: [(name, xyz)]}) in residue order."""
    ox = _cage_oxygens_by_residue(seq)
    serial = 1
    with open(path, "w") as fh:
        for i, aa in enumerate(seq):
            t = i * 0.7
            ca = (20.0 * np.cos(t), 20.0 * np.sin(t), 0.0)
            # extra_atoms may override the CA position for a planted residue.
            planted = {n: xyz for n, xyz in extra_atoms.get(i, [])}
            if "CA" not in planted:
                fh.write(_atom_line(serial, "CA", _THREE[aa], "A", i + 1, ca)); serial += 1
            for name, xyz in extra_atoms.get(i, []):
                fh.write(_atom_line(serial, name, _THREE[aa], "A", i + 1, xyz)); serial += 1
            for name, xyz in ox.get(i, []):
                fh.write(_atom_line(serial, name, _THREE[aa], "A", i + 1, xyz)); serial += 1
        fh.write("END\n")


def test_counts_and_RY_pair_by_adjacency():
    # Arg(idx0): NH1 near origin (d=2), CA farther (d=5) -> reaches toward site -> counts.
    # Lys(idx1): NZ near (d=3), CA farther (d=8) -> counts.
    # Tyr(idx2): within +-2 residues of the Arg -> RY pair by sequence adjacency.
    extra = {
        0: [("CA", (5.0, 0, 0)), ("CZ", (2.2, 0, 0)), ("NH1", (2.0, 0, 0))],
        1: [("CA", (0, 8.0, 0)), ("NZ", (0, 3.0, 0))],
        2: [("CA", (0, 0, 9.0)), ("OH", (30.0, 0, 0))],  # OH far: adjacency, not spatial
    }
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write(p, _SEQ, extra)
        r = diphosphate_sensor_one(ResidueInfo(p))
        assert r["metal_point_found"] is True
        assert r["n_arg"] == 1, r
        assert r["n_lys"] == 1, r
        assert r["n_diphosphate_basic_residues"] == 2
        assert r["has_RY_pair"] is True
        assert r["n_RY_pairs"] == 1
        assert r["n_residues"] == len(_SEQ)


def test_RY_pair_by_spatial_proximity():
    # Tyr is NOT adjacent to the Arg in sequence, but its OH is close to the Arg
    # guanidinium centroid AND near the metal point -> spatial RY pair.
    extra = {
        0: [("CA", (5.0, 0, 0)), ("CZ", (2.2, 0, 0)), ("NH1", (2.0, 0, 0))],
        20: [("CA", (0, 0, 9.0)), ("OH", (2.5, 0, 0))],   # near arg centroid (~2.0,0,0) & site
    }
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write(p, _SEQ, extra)
        r = diphosphate_sensor_one(ResidueInfo(p))
        assert r["n_arg"] == 1
        assert r["has_RY_pair"] is True
        assert r["n_RY_pairs"] == 1


def test_direction_and_cutoff_exclusions():
    # Arg(idx0): terminal N far (d=10) but CA nearer (d=1) -> side chain points AWAY -> excluded.
    # Lys(idx1): NZ beyond the 12 A cutoff -> excluded.
    extra = {
        0: [("CA", (1.0, 0, 0)), ("NH1", (10.0, 0, 0))],
        1: [("CA", (0, 25.0, 0)), ("NZ", (0, 20.0, 0))],
    }
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write(p, _SEQ, extra)
        r = diphosphate_sensor_one(ResidueInfo(p))
        assert r["metal_point_found"] is True
        assert r["n_arg"] == 0, r
        assert r["n_lys"] == 0, r
        assert r["n_diphosphate_basic_residues"] == 0
        assert r["has_RY_pair"] is False


def test_no_motif_no_metal_point():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "poly.pdb")
        with open(p, "w") as fh:
            for i in range(30):
                fh.write(_atom_line(i + 1, "CA", "ALA", "A", i + 1, (i * 3.8, 0, 0)))
            fh.write("END\n")
        r = diphosphate_sensor_one(ResidueInfo(p))
        assert r["metal_point_found"] is False
        assert r["n_diphosphate_basic_residues"] == 0
        assert r["has_RY_pair"] is False
        assert r["n_residues"] == 30


def test_dir_id_keying_csv_and_nan_on_broken():
    extra = {
        0: [("CA", (5.0, 0, 0)), ("CZ", (2.2, 0, 0)), ("NH1", (2.0, 0, 0))],
        1: [("CA", (0, 8.0, 0)), ("NZ", (0, 3.0, 0))],
        2: [("CA", (0, 0, 9.0)), ("OH", (30.0, 0, 0))],
    }
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write(os.path.join(structs, "good.pdb"), _SEQ, extra)
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("not a pdb\n")
        df = diphosphate_sensor_dir(structs)

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        assert os.path.isfile(structs + "_diphosphate_sensor.csv")

        good = df.set_index("ID").loc["good"]
        assert bool(good["metal_point_found"]) is True
        assert int(good["n_diphosphate_basic_residues"]) == 2

        broken = df.set_index("ID").loc["broken"]
        assert bool(broken["metal_point_found"]) is False
        assert int(broken["n_residues"]) == 0


def test_default_cutoff_value():
    # Guard the documented default so a silent change is caught.
    assert DEFAULT_CUTOFF == 12.0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
