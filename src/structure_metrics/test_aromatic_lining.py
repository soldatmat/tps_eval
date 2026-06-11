from __future__ import annotations

"""Self-contained tests for aromatic_lining.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_aromatic_lining.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_aromatic_lining.py -q

Strategy: build tiny synthetic PDBs whose SEQUENCE contains a DDXXD motif and an
NSE/DTE motif (so the shared localizer + the reused carboxylate-cage machinery can
place a metal point at a known location), then plant aromatic residues at controlled
distances and ring orientations relative to the locus. The geometry is chosen so the
expected counts are known exactly. Also exercises the no-metal-point red-flag row and
the per-directory ID keying / NaN contract.
"""

import os
import tempfile

import numpy as np

from aromatic_lining import (
    COLUMNS,
    _ring_centroid_normal,
    aromatic_lining,
    aromatic_lining_dir,
)


# --- minimal-PDB construction helpers ---------------------------------------

# Map 1-letter -> 3-letter for the residue types we emit.
_THREE = {
    "A": "ALA", "D": "ASP", "E": "GLU", "N": "ASN", "S": "SER", "T": "THR",
    "L": "LEU", "I": "ILE", "V": "VAL", "G": "GLY",
    "F": "PHE", "Y": "TYR", "W": "TRP", "H": "HIS",
}


class _PDBWriter:
    def __init__(self):
        self.lines = []
        self.serial = 1

    def atom(self, resname, resseq, atomname, xyz, element=None):
        x, y, z = xyz
        el = element or atomname[0]
        # columns per PDB spec; atom name left-justified in cols 13-16 when 4 chars,
        # else padded with a leading space (cols 14-16). Biopython tolerates the
        # simple " {name:<3}" form for our 1-2 char names.
        # PDB cols: 13-16 name, 17 altLoc, 18-20 resName, 21 blank, 22 chainID,
        # 23-26 resSeq. Atom names <4 chars are right-justified from col 14.
        name_field = f" {atomname:<3}"[:4]
        self.lines.append(
            f"ATOM  {self.serial:>5d} {name_field} {resname:>3} A{resseq:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}"
        )
        self.serial += 1

    def write(self, path):
        with open(path, "w") as fh:
            fh.write("\n".join(self.lines) + "\nEND\n")


def _planar_ring(centroid, normal, radius=1.39):
    """6 coords of a regular hexagon centred at ``centroid`` in the plane whose
    normal is ``normal``. Returns coords for (CG/CD2 .. ) in ring order."""
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    # two in-plane orthonormal vectors
    a = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(a, normal)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    u = a - np.dot(a, normal) * normal
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    pts = []
    for k in range(6):
        ang = 2 * np.pi * k / 6
        pts.append(np.asarray(centroid) + radius * (np.cos(ang) * u + np.sin(ang) * v))
    return pts


def _add_backbone(w, resname, resseq, ca_xyz):
    """Add a minimal backbone (just CA, plus CB for non-Gly) at ca_xyz."""
    w.atom(resname, resseq, "CA", ca_xyz, "C")
    if resname != "GLY":
        w.atom(resname, resseq, "CB", np.asarray(ca_xyz) + np.array([0.5, 0.0, 0.0]), "C")


def _add_coordinating_oxygens(w, resname, resseq, ca_xyz, o_target):
    """Place the residue's side-chain coordinating oxygen(s) AT (near) ``o_target``
    so the carboxylate-cage centroid lands where we want."""
    t = np.asarray(o_target, dtype=float)
    if resname == "ASP":
        w.atom(resname, resseq, "OD1", t + np.array([0.0, 0.05, 0.0]), "O")
        w.atom(resname, resseq, "OD2", t - np.array([0.0, 0.05, 0.0]), "O")
    elif resname == "GLU":
        w.atom(resname, resseq, "OE1", t + np.array([0.0, 0.05, 0.0]), "O")
        w.atom(resname, resseq, "OE2", t - np.array([0.0, 0.05, 0.0]), "O")
    elif resname == "ASN":
        w.atom(resname, resseq, "OD1", t, "O")
    elif resname == "SER":
        w.atom(resname, resseq, "OG", t, "O")
    elif resname == "THR":
        w.atom(resname, resseq, "OG1", t, "O")


def _add_aromatic_ring(w, resname, resseq, centroid, normal):
    """Place the six-membered ring atoms (names per RING_ATOMS) for an aromatic."""
    ring = _planar_ring(centroid, normal)
    if resname in ("PHE", "TYR"):
        names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
    else:  # TRP
        names = ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")
    for nm, xyz in zip(names, ring):
        w.atom(resname, resseq, nm, xyz, "C")


def _build_motif_structure(path, metal_point, aromatics):
    """Build a PDB whose derived sequence has a DDXXD + NSE/DTE motif with their
    coordinating oxygens clustered at ``metal_point``, plus the given aromatics.

    ``aromatics`` is a list of dicts: {resname, ca, ring_centroid, ring_normal}.
    Sequence layout (positions 1-based):
      1-5   DDAAD       (DDXXD: coords at 1,2,5)
      6-10  AAAAA       spacer
      11-19 NDLASAAAE   (NSE/DTE (N|D)D(L|I|V).(S|T)...E: coords at 11,12,15,19)
      then  one residue per aromatic, then 'A' padding.

    The motif/spacer CAs are scattered SYMMETRICALLY on a sphere around the metal
    point so the pocket-residue Calpha centroid (the cation-locus approximation)
    lands ~at the metal point, keeping the test geometry deterministic.
    """
    w = _PDBWriter()
    mp = np.asarray(metal_point, dtype=float)

    # 19 symmetric directions on a sphere (Fibonacci) for the 5+5+9 motif/spacer CAs.
    def _sphere(n, radius=8.0):
        pts = []
        for k in range(n):
            z = 1.0 - 2.0 * (k + 0.5) / n
            r = np.sqrt(max(0.0, 1.0 - z * z))
            phi = np.pi * (3.0 - np.sqrt(5.0)) * k
            pts.append(mp + radius * np.array([r * np.cos(phi), r * np.sin(phi), z]))
        return pts

    dirs = _sphere(19)

    # DDXXD: D D A A D  -> coordinating at offsets 0,1,4
    seq1 = "DDAAD"
    for i, aa in enumerate(seq1):
        resseq = i + 1
        rn = _THREE[aa]
        ca = dirs[i]
        _add_backbone(w, rn, resseq, ca)
        if aa in ("D",):
            _add_coordinating_oxygens(w, rn, resseq, ca, mp)

    # spacer AAAAA
    for i in range(5):
        resseq = 6 + i
        _add_backbone(w, "ALA", resseq, dirs[5 + i])

    # NSE/DTE: N D L A S A A A E -> coordinating at offsets 0,1,4,8
    seq2 = "NDLASAAAE"
    coord_off = {0, 1, 4, 8}
    for i, aa in enumerate(seq2):
        resseq = 11 + i
        rn = _THREE[aa]
        ca = dirs[10 + i]
        _add_backbone(w, rn, resseq, ca)
        if i in coord_off:
            _add_coordinating_oxygens(w, rn, resseq, ca, mp)

    # aromatics
    next_seq = 11 + len(seq2)
    for arom in aromatics:
        rn = arom["resname"]
        ca = np.asarray(arom["ca"], dtype=float)
        _add_backbone(w, rn, next_seq, ca)
        _add_aromatic_ring(w, rn, next_seq, arom["ring_centroid"], arom["ring_normal"])
        next_seq += 1

    w.write(path)


# --- tests -------------------------------------------------------------------


def test_ring_centroid_normal_recovers_plane():
    """The SVD normal of a planar hexagon equals the plane normal (up to sign)."""
    from Bio.PDB import PDBParser
    import warnings
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    centroid = np.array([10.0, 0.0, 0.0])
    normal = np.array([0.0, 0.0, 1.0])
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "phe.pdb")
        w = _PDBWriter()
        _add_backbone(w, "PHE", 1, centroid + np.array([0, 0, -2]))
        _add_aromatic_ring(w, "PHE", 1, centroid, normal)
        w.write(p)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = PDBParser(QUIET=True).get_structure("s", p)
        res = next(iter(next(iter(next(iter(s))))))
        got_c, got_n = _ring_centroid_normal(res)
    np.testing.assert_allclose(got_c, centroid, atol=1e-2)
    assert abs(abs(float(np.dot(got_n, normal))) - 1.0) < 1e-3


def test_counts_and_inward_facing():
    """Two pocket aromatics (a Trp + a Phe), one far Tyr (outside cutoff). The Trp is
    placed face-on within cation-pi range; the Phe is edge-on (its normal is
    perpendicular to the locus direction) so it should NOT count as inward."""
    metal_point = np.array([0.0, 0.0, 0.0])
    # locus = CA centroid of pocket residues. The motif residues' CAs surround the
    # origin, so the locus is near origin; we keep the aromatic ring centroids close
    # so distance to the locus stays in the 3.5-6 A window.
    aromatics = [
        # In-pocket TRP, ring centroid ~4 A out along +x, FACE pointing back toward
        # the locus (normal ~ along x) -> inward.
        {
            "resname": "TRP", "ca": np.array([6.0, 0.0, 0.0]),
            "ring_centroid": np.array([4.0, 0.0, 0.0]),
            "ring_normal": np.array([1.0, 0.0, 0.0]),
        },
        # In-pocket PHE, ring centroid ~4 A out along +y, EDGE-on (normal along z,
        # perpendicular to the y direction toward locus) -> NOT inward.
        {
            "resname": "PHE", "ca": np.array([0.0, 6.0, 0.0]),
            "ring_centroid": np.array([0.0, 4.0, 0.0]),
            "ring_normal": np.array([0.0, 0.0, 1.0]),
        },
        # FAR TYR well outside the 10 A cutoff -> not a pocket residue at all.
        {
            "resname": "TYR", "ca": np.array([40.0, 40.0, 40.0]),
            "ring_centroid": np.array([40.0, 40.0, 42.0]),
            "ring_normal": np.array([0.0, 0.0, 1.0]),
        },
    ]
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "design.pdb")
        _build_motif_structure(p, metal_point, aromatics)
        r = aromatic_lining(p, cutoff=10.0)

    assert r["metal_point_found"] is True
    assert r["n_trp"] == 1, r
    assert r["n_phe"] == 1, r
    assert r["n_tyr"] == 0, r  # far Tyr excluded by cutoff
    assert r["n_pocket_aromatics"] == 2, r
    assert 0.0 < r["aromatic_fraction"] <= 1.0
    # Only the face-on Trp counts as inward; the edge-on Phe does not.
    assert r["n_inward_facing_aromatics"] == 1, r


def test_no_metal_point_is_red_flag():
    """A structure whose sequence lacks the motifs gets metal_point_found False and
    NaN counts (not dropped)."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "nomotif.pdb")
        w = _PDBWriter()
        for i in range(10):
            _add_backbone(w, "ALA", i + 1, np.array([float(i), 0.0, 0.0]))
        w.write(p)
        r = aromatic_lining(p)
    assert r["metal_point_found"] is False
    assert np.isnan(r["n_pocket_aromatics"])
    assert np.isnan(r["aromatic_fraction"])
    assert r["n_residues"] == 10


def test_dir_keyed_by_id_and_columns():
    metal_point = np.array([0.0, 0.0, 0.0])
    aromatics = [
        {
            "resname": "TRP", "ca": np.array([6.0, 0.0, 0.0]),
            "ring_centroid": np.array([4.0, 0.0, 0.0]),
            "ring_normal": np.array([1.0, 0.0, 0.0]),
        },
    ]
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _build_motif_structure(os.path.join(structs, "good.pdb"), metal_point, aromatics)
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("not a pdb\n")
        df = aromatic_lining_dir(structs)
        # The sibling CSV is written next to the dir (check before the temp dir is removed).
        assert os.path.isfile(structs + "_aromatic_lining.csv")

    assert list(df.columns) == COLUMNS
    assert set(df["ID"]) == {"good", "broken"}
    good = df.set_index("ID").loc["good"]
    assert bool(good["metal_point_found"]) is True
    assert int(good["n_trp"]) == 1
    broken = df.set_index("ID").loc["broken"]
    # broken parses to an empty structure -> no motif -> red flag.
    assert bool(broken["metal_point_found"]) is False


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
