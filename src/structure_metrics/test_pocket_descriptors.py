from __future__ import annotations

"""Self-contained tests for pocket_descriptors.py (numpy + pandas + biopython only).

Run from this directory so the flat-module imports resolve:
    cd src/structure_metrics && python test_pocket_descriptors.py
or:
    cd src/structure_metrics && python -m pytest test_pocket_descriptors.py -q

fpocket and P2Rank are EXTERNAL binaries not present in this env, so the engine calls
degrade gracefully to NaN — which is itself part of the contract and is tested here.
The engine-independent logic IS exercised on synthetic inputs: the metal-point anchor
(closed-form origin from a two-motif carboxylate cage), the fpocket ``*_info.txt`` and
alpha-sphere ``*.pqr`` parsers (incl. the "Volume" vs "Volume score" disambiguation),
the P2Rank predictions-CSV parser + rank assignment, the bounding-box enclosure test,
the ``pocket_sasa_per_volume`` derivation (via a monkeypatched fpocket result), the
metal-point-absent -> all-NaN contract, ID keying, the sibling-CSV filename, and
NaN-on-broken-structure.
"""

import os
import tempfile

import numpy as np

import pocket_descriptors as pk
from pocket_descriptors import (
    COLUMNS,
    METAL_POINT_CUTOFF_A,
    _nan_result,
    _parse_fpocket_alpha_spheres,
    _parse_fpocket_info,
    _parse_p2rank_predictions,
    _point_inside_cloud,
    metal_point,
    pocket_descriptors,
    pocket_descriptors_dir,
)

_THREE = {"A": "ALA", "D": "ASP", "N": "ASN", "L": "LEU", "S": "SER", "E": "GLU"}
_SEQ = "AAAAA" + "DDAAD" + "AAAAA" + "NDLASAAAE" + "AA"
_COORD_IDX = [5, 6, 9, 15, 16, 19, 23]
_OX_ATOMS = {"ASP": ("OD1", "OD2"), "ASN": ("OD1",), "SER": ("OG",), "GLU": ("OE1", "OE2")}
_R = 3.0
_AXIS6 = [(_R, 0, 0), (-_R, 0, 0), (0, _R, 0), (0, -_R, 0), (0, 0, _R), (0, 0, -_R)]
_SPHERE12 = _AXIS6 + _AXIS6


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _atom_line(serial, name, resname, chain, resseq, xyz, *, element=None):
    x, y, z = xyz
    if element is None:
        element = "O" if name.startswith("O") else name[0]
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{'ATOM':<6}{serial:>5} {atom_field}{' ':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}{'':10}{element:>2}\n"
    )


def _write_cage_pdb(path, seq=_SEQ):
    """Two-motif carboxylate cage: 12 coordinating oxygens summing to zero -> metal
    point at the origin. Every residue gets a far-away CA (radius 20)."""
    slots = []
    for i in _COORD_IDX:
        for name in _OX_ATOMS[_THREE[seq[i]]]:
            slots.append((i, name))
    ox = {}
    for (i, name), xyz in zip(slots, _SPHERE12):
        ox.setdefault(i, []).append((name, xyz))
    serial = 1
    with open(path, "w") as fh:
        for i, aa in enumerate(seq):
            t = i * 0.7
            fh.write(_atom_line(serial, "CA", _THREE[aa], "A", i + 1,
                                (20.0 * np.cos(t), 20.0 * np.sin(t), 0.0), element="C"))
            serial += 1
            for name, xyz in ox.get(i, []):
                fh.write(_atom_line(serial, name, _THREE[aa], "A", i + 1, xyz)); serial += 1
        fh.write("END\n")


# ------------------------------ metal point ----------------------------------- #
def test_metal_point_origin_from_cage():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "cage.pdb")
        _write_cage_pdb(p)
        point, n = metal_point(p)
        assert point is not None
        np.testing.assert_allclose(point, [0.0, 0.0, 0.0], atol=1e-6)
        assert n == len(_SEQ)


def test_metal_point_none_without_motif():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "poly.pdb")
        with open(p, "w") as fh:
            for i in range(30):
                fh.write(_atom_line(i + 1, "CA", "ALA", "A", i + 1, (i * 3.8, 0, 0), element="C"))
            fh.write("END\n")
        point, n = metal_point(p)
        assert point is None
        assert n == 30


# ------------------------------ fpocket parsers ------------------------------- #
def test_parse_fpocket_alpha_spheres():
    with tempfile.TemporaryDirectory() as d:
        pqr = os.path.join(d, "pocket1_vert.pqr")
        with open(pqr, "w") as fh:
            fh.write(_atom_line(1, "APOL", "STP", "A", 1, (1, 2, 3), element="C"))
            fh.write(_atom_line(2, "APOL", "STP", "A", 1, (4, 5, 6), element="C"))
            fh.write(_atom_line(3, "CA", "ALA", "A", 1, (9, 9, 9), element="C"))  # non-STP ignored
        cloud = _parse_fpocket_alpha_spheres(pqr)
        assert cloud.shape == (2, 3)
        np.testing.assert_allclose(cloud, [[1, 2, 3], [4, 5, 6]])
    # Missing file -> empty, not a crash.
    assert _parse_fpocket_alpha_spheres("/no/such.pqr").shape == (0, 3)


def test_parse_fpocket_info_volume_vs_volume_score():
    info = """Pocket 1 :
\tScore : \t\t0.5
\tVolume : \t\t1000.0
\tVolume score : \t\t4.5
\tNumber of Alpha Spheres : \t50
\tTotal SASA : \t200.0
\tHydrophobicity score : \t30.0
\tMean local hydrophobic density : \t12.0

Pocket 2 :
\tVolume : \t\t500.0
"""
    with tempfile.TemporaryDirectory() as d:
        info_path = os.path.join(d, "x_info.txt")
        with open(info_path, "w") as fh:
            fh.write(info)
        pockets = _parse_fpocket_info(info_path)
        assert set(pockets) == {1, 2}
        p1 = pockets[1]
        # "Volume score" must NOT overwrite the exact "Volume".
        _approx(p1["catalytic_pocket_volume"], 1000.0)
        _approx(p1["pocket_n_alpha_spheres"], 50.0)
        _approx(p1["pocket_total_sasa"], 200.0)
        _approx(p1["pocket_hydrophobicity"], 30.0)
        _approx(p1["pocket_enclosure"], 12.0)
        _approx(pockets[2]["catalytic_pocket_volume"], 500.0)
    assert _parse_fpocket_info("/no/such_info.txt") == {}


# ------------------------------ P2Rank parser --------------------------------- #
def test_parse_p2rank_predictions_rank_by_file_order():
    csv = ("name, rank, score, center_x, center_y, center_z\n"
           "pocket1, 1, 0.90, 1.0, 0.0, 0.0\n"
           "pocket2, 2, 0.50, 10.0, 0.0, 0.0\n")
    with tempfile.TemporaryDirectory() as d:
        pred = os.path.join(d, "x.pdb_predictions.csv")
        with open(pred, "w") as fh:
            fh.write(csv)
        pockets = _parse_p2rank_predictions(pred)
        assert len(pockets) == 2
        assert pockets[0]["rank"] == 1 and pockets[1]["rank"] == 2
        _approx(pockets[0]["score"], 0.90)
        np.testing.assert_allclose(pockets[0]["center"], [1.0, 0.0, 0.0])
    # A CSV missing required columns -> empty list (no crash).
    with tempfile.TemporaryDirectory() as d:
        bad = os.path.join(d, "bad.csv")
        with open(bad, "w") as fh:
            fh.write("a,b\n1,2\n")
        assert _parse_p2rank_predictions(bad) == []


# ------------------------------ enclosure test -------------------------------- #
def test_point_inside_cloud():
    cloud = np.array([[0, 0, 0], [10, 10, 10]], dtype=float)
    assert _point_inside_cloud(np.array([5.0, 5, 5]), cloud) is True
    assert _point_inside_cloud(np.array([50.0, 5, 5]), cloud) is False
    # Padding admits a point just outside the raw bounding box.
    assert _point_inside_cloud(np.array([-2.0, 0, 0]), cloud, pad=3.0) is True
    assert _point_inside_cloud(np.array([0.0, 0, 0]), np.empty((0, 3))) is False


# ------------------------------ _nan_result ----------------------------------- #
def test_nan_result_shape():
    r = _nan_result(17)
    assert r["n_residues"] == 17
    assert r["metal_point_found"] is False
    assert r["fpocket_catalytic_pocket_found"] is False
    assert r["p2rank_catalytic_pocket_found"] is False
    for k in ("catalytic_pocket_volume", "pocket_total_sasa", "pocket_sasa_per_volume",
              "p2rank_catalytic_site_score"):
        assert np.isnan(r[k]), k


# ---------------------- graceful degradation (no fpocket) --------------------- #
def test_pocket_descriptors_no_fpocket_binary():
    # fpocket is absent -> run_fpocket returns None -> all fpocket columns NaN, but the
    # metal point is still found and the row is well-formed (no crash).
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "cage.pdb")
        _write_cage_pdb(p)
        r = pocket_descriptors(p, fpocket_bin="fpocket_definitely_missing_bin", p2rank_bin=None)
        assert r["metal_point_found"] is True
        assert r["fpocket_catalytic_pocket_found"] is False
        assert np.isnan(r["catalytic_pocket_volume"])
        assert np.isnan(r["pocket_sasa_per_volume"])
        assert np.isnan(r["p2rank_catalytic_site_score"])   # p2rank skipped
        assert r["n_residues"] == len(_SEQ)


def test_pocket_descriptors_no_motif_all_nan():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "poly.pdb")
        with open(p, "w") as fh:
            for i in range(20):
                fh.write(_atom_line(i + 1, "CA", "ALA", "A", i + 1, (i * 3.8, 0, 0), element="C"))
            fh.write("END\n")
        r = pocket_descriptors(p, p2rank_bin=None)
        assert r["metal_point_found"] is False
        assert np.isnan(r["catalytic_pocket_volume"])


def test_sasa_per_volume_derivation_via_monkeypatch():
    # Exercise the derived specific-surface-area formula (Total SASA / Volume) by
    # substituting a known fpocket result. Restore the original after.
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "cage.pdb")
        _write_cage_pdb(p)
        orig = pk.fpocket_catalytic
        pk.fpocket_catalytic = lambda *a, **k: {
            "catalytic_pocket_volume": 1000.0,
            "pocket_hydrophobicity": np.nan,
            "pocket_enclosure": np.nan,
            "pocket_n_alpha_spheres": np.nan,
            "pocket_total_sasa": 250.0,
            "pocket_depth": np.nan,
            "fpocket_catalytic_pocket_found": True,
        }
        try:
            r = pocket_descriptors(p, p2rank_bin=None)
        finally:
            pk.fpocket_catalytic = orig
        _approx(r["pocket_sasa_per_volume"], 0.25)
        # Volume <= 0 -> NaN, not a division error.
        pk.fpocket_catalytic = lambda *a, **k: {
            "catalytic_pocket_volume": 0.0, "pocket_total_sasa": 250.0,
            "fpocket_catalytic_pocket_found": True,
        }
        try:
            r0 = pocket_descriptors(p, p2rank_bin=None)
        finally:
            pk.fpocket_catalytic = orig
        assert np.isnan(r0["pocket_sasa_per_volume"])


# ------------------------------ dir driver ------------------------------------ #
def test_dir_id_keying_csv_and_nan_on_broken():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write_cage_pdb(os.path.join(structs, "good.pdb"))
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("not a pdb\n")
        df = pocket_descriptors_dir(structs, fpocket_bin="fpocket_missing_bin")

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        assert os.path.isfile(structs + "_pocket_descriptors.csv")

        good = df.set_index("ID").loc["good"]
        assert bool(good["metal_point_found"]) is True    # motif found even w/o fpocket
        broken = df.set_index("ID").loc["broken"]
        assert bool(broken["metal_point_found"]) is False
        assert int(broken["n_residues"]) == 0


def test_cutoff_constant():
    assert METAL_POINT_CUTOFF_A == 12.0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
