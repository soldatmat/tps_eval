from __future__ import annotations

"""Self-contained tests for motif_structural_distance.py (numpy + biopython only — no
conda env / EnzymeExplorer / PyMOL). Run from this directory so the flat-module imports
resolve:
    cd src/structure_metrics && python test_motif_structural_distance.py
or:
    cd src/structure_metrics && python -m pytest test_motif_structural_distance.py -q

Builds a synthetic protein whose residue NAMES spell a sequence carrying exactly one
DDXXD-family motif and one NSE/DTE motif at known positions, with the metal-coordinating
CA atoms placed at closed-form coordinates so ``motif_centroid_distance`` (== 10 A) and
``motif_min_ca_distance`` (== 8 A) are exact. Also covers the missing-motif -> NaN
contract, HETATM/missing-CA handling in the sequence/CA extraction, ID keying, the
sibling-CSV filename, and NaN-on-broken-structure.
"""

import os
import tempfile

import numpy as np

from motif_structural_distance import (
    COLUMNS,
    _coord_matrix,
    motif_distances,
    motif_structural_distance_dir,
    structure_sequence_and_ca,
)

# 1-letter -> 3-letter for the residue names we use.
_THREE = {"A": "ALA", "D": "ASP", "N": "ASN", "L": "LEU", "S": "SER", "E": "GLU",
          "G": "GLY"}


def _approx(a, b, tol=1e-3):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _atom_line(serial, name, resname, chain, resseq, xyz, *, record="ATOM", element="C"):
    x, y, z = xyz
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{record:<6}{serial:>5} {atom_field}{' ':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}{'':10}{element:>2}\n"
    )


def _write_residues(path, residues, *, chain="A"):
    """residues: list of (one_letter, ca_xyz | None). One CA per residue (skipped when
    ca_xyz is None, so we can exercise the missing-CA branch)."""
    serial = 1
    with open(path, "w") as fh:
        for i, (aa, xyz) in enumerate(residues, start=1):
            if xyz is not None:
                fh.write(_atom_line(serial, "CA", _THREE[aa], chain, i, xyz))
                serial += 1
        fh.write("END\n")


# Sequence with the two motifs at known positions (verified against motif_localization):
#   DDXXD 'DDAAD' at index 5 (coordinating 5,6,9),
#   NSE/DTE 'NDLASAAAE' at index 15 (coordinating 15,16,19,23).
_SEQ = "AAAAA" + "DDAAD" + "AAAAA" + "NDLASAAAE" + "AA"

# CA coordinates: coordinating CAs placed so the DDXXD centroid is the origin and the
# NSE/DTE centroid is (10,0,0); every other residue is parked far away.
_COORDS = {5: (-1, 0, 0), 6: (1, 0, 0), 9: (0, 0, 0),
           15: (10, -1, 0), 16: (10, 1, 0), 19: (11, 0, 0), 23: (9, 0, 0)}


def _motif_residues():
    return [(_SEQ[i], _COORDS.get(i, (0.0, 0.0, 50.0 + i))) for i in range(len(_SEQ))]


def test_known_centroid_and_min_distance():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write_residues(p, _motif_residues())
        r = motif_distances(p)
        assert r["n_residues"] == len(_SEQ)
        _approx(r["motif_centroid_distance"], 10.0)
        _approx(r["motif_min_ca_distance"], 8.0)


def test_missing_motif_is_nan():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "poly_ala.pdb")
        _write_residues(p, [("A", (i * 3.8, 0, 0)) for i in range(30)])
        r = motif_distances(p)
        assert np.isnan(r["motif_centroid_distance"])
        assert np.isnan(r["motif_min_ca_distance"])
        assert r["n_residues"] == 30


def test_ddxxd_only_is_nan():
    # DDXXD present but no NSE/DTE -> both distances NaN (needs BOTH motifs).
    seq = "AAAAA" + "DDAAD" + "AAAAAAAAAA"
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "one_motif.pdb")
        _write_residues(p, [(aa, (i * 3.8, 0, 0)) for i, aa in enumerate(seq)])
        r = motif_distances(p)
        assert np.isnan(r["motif_centroid_distance"])


def test_structure_sequence_and_ca_skips_hetatm_and_missing_ca():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "s.pdb")
        with open(p, "w") as fh:
            fh.write(_atom_line(1, "CA", "ALA", "A", 1, (0, 0, 0)))
            fh.write(_atom_line(2, "N", "GLY", "A", 2, (4, 0, 0), element="N"))  # residue w/o CA
            fh.write(_atom_line(3, "CA", "ASP", "A", 3, (8, 0, 0)))
            # HETATM ligand carbon must be skipped (not part of the sequence).
            fh.write(_atom_line(999, "C1", "FPP", "A", 900, (99, 99, 99),
                                record="HETATM", element="C"))
            fh.write("END\n")
        seq, ca = structure_sequence_and_ca(p)
        assert seq == "AGD", seq                 # HETATM excluded
        assert ca[0] is not None and ca[2] is not None
        assert ca[1] is None                     # GLY had no CA -> None, index-aligned


def test_nonstandard_residue_maps_to_X():
    # A non-standard ATOM residue (e.g. UNK) must degrade to 'X', not crash the parse.
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "unk.pdb")
        with open(p, "w") as fh:
            fh.write(_atom_line(1, "CA", "ALA", "A", 1, (0, 0, 0)))
            fh.write(_atom_line(2, "CA", "UNK", "A", 2, (4, 0, 0)))
            fh.write("END\n")
        seq, _ = structure_sequence_and_ca(p)
        assert seq == "AX", seq


def test_coord_matrix_drops_missing_and_out_of_range():
    ca = [np.array([0.0, 0, 0]), None, np.array([3.0, 0, 0])]
    m = _coord_matrix([0, 1, 2, 99], ca)   # index 1 is None, 99 out of range
    assert m.shape == (2, 3)
    assert _coord_matrix([1], ca) is None   # only a missing CA -> None


def test_dir_id_keying_csv_and_nan_on_broken():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write_residues(os.path.join(structs, "good.pdb"), _motif_residues())
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("not a pdb\n")
        df = motif_structural_distance_dir(structs)

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        assert os.path.isfile(structs + "_motif_structural_distance.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["motif_centroid_distance"]), 10.0)
        broken = df.set_index("ID").loc["broken"]
        assert np.isnan(broken["motif_centroid_distance"])
        assert int(broken["n_residues"]) == 0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
