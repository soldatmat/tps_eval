from __future__ import annotations

"""Self-contained tests for active_site_geometry.py (numpy + biopython only — no conda
env / EnzymeExplorer / PyMOL). Run from this directory so the flat-module imports
resolve:
    cd src/structure_metrics && python test_active_site_geometry.py
or:
    cd src/structure_metrics && python -m pytest test_active_site_geometry.py -q

Builds a synthetic protein whose residue names spell one DDXXD + one NSE/DTE motif, with
the metal-coordinating side-chain oxygens placed at 12 points forming three antipodal
axis pairs at radius 3 about the origin. That makes the geometry closed-form:
``carboxylate_convergence_radius`` == 3, oxygen centroid (metal point) == origin, and
with every CA parked on a radius-5 sphere the ``metal_point_void`` == 5. Also covers the
relaxed DDXXD-only path, the missing-motif -> NaN contract, the constellation-RMSD
template path (translation-invariant ~0 on a rigidly shifted copy), ID keying, the
sibling-CSV filename, and NaN-on-broken-structure.
"""

import os
import tempfile

import numpy as np

from active_site_geometry import (
    COLUMNS,
    COORDINATING_OXYGEN_ATOMS,
    _coordinating_indices_both,
    active_site_geometry,
    active_site_geometry_dir,
    build_templates,
    coordinating_indices_relaxed,
    metal_point,
    structure_sequence_residues_atoms,
)

_THREE = {"A": "ALA", "D": "ASP", "N": "ASN", "L": "LEU", "S": "SER", "E": "GLU"}

# Same two-motif layout used by test_motif_structural_distance (verified positions):
#   DDXXD 'DDAAD' at 5 (coordinating 5,6,9), NSE/DTE 'NDLASAAAE' at 15 (16,19,23 too).
_SEQ = "AAAAA" + "DDAAD" + "AAAAA" + "NDLASAAAE" + "AA"
_COORD_IDX = [5, 6, 9, 15, 16, 19, 23]

# 12 oxygen positions = three antipodal axis pairs at r=3, each used twice -> they sum
# to zero (centroid == origin) and each sits at distance 3 (RMS radius == 3).
_R = 3.0
_AXIS6 = [(_R, 0, 0), (-_R, 0, 0), (0, _R, 0), (0, -_R, 0), (0, 0, _R), (0, 0, -_R)]
_SPHERE12 = _AXIS6 + _AXIS6
_CA_RADIUS = 5.0


def _approx(a, b, tol=1e-3):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _atom_line(serial, name, resname, chain, resseq, xyz, *, element="C"):
    x, y, z = xyz
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{'ATOM':<6}{serial:>5} {atom_field}{' ':1}{resname:>3} "
        f"{chain}{resseq:>4}{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}"
        f"{'':10}{element:>2}\n"
    )


def _cage_atoms(seq=_SEQ, *, shift=(0.0, 0.0, 0.0)):
    """Return the atom list (name, resname, resseq, xyz) for the two-motif cage. Every
    residue gets a CA on a radius-5 circle; the coordinating residues additionally get
    their side-chain oxygens at the radius-3 sphere points. ``shift`` rigidly translates
    everything (for translation-invariance tests)."""
    sx, sy, sz = shift
    atoms = []
    # Enumerate the coordinating (residue_index, oxygen_atom_name) slots in gather order.
    slots = []
    for i in _COORD_IDX:
        for name in COORDINATING_OXYGEN_ATOMS[_THREE[seq[i]]]:
            slots.append((i, name))
    assert len(slots) == len(_SPHERE12), (len(slots), len(_SPHERE12))
    ox_by_res = {}
    for (i, name), (ox, oy, oz) in zip(slots, _SPHERE12):
        ox_by_res.setdefault(i, []).append((name, (ox + sx, oy + sy, oz + sz)))

    for i, aa in enumerate(seq):
        t = i * 0.5
        ca = (_CA_RADIUS * np.cos(t) + sx, _CA_RADIUS * np.sin(t) + sy, sz)
        atoms.append(("CA", _THREE[aa], i + 1, ca))
        for name, xyz in ox_by_res.get(i, []):
            atoms.append((name, _THREE[aa], i + 1, xyz))
    return atoms


def _write_atoms(path, atoms, *, chain="A"):
    with open(path, "w") as fh:
        for serial, (name, resname, resseq, xyz) in enumerate(atoms, start=1):
            element = "O" if name.startswith("O") else name[0]
            fh.write(_atom_line(serial, name, resname, chain, resseq, xyz, element=element))
        fh.write("END\n")


def _write_poly_ala(path, n=30):
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(_atom_line(i + 1, "CA", "ALA", "A", i + 1, (i * 3.8, 0, 0)))
        fh.write("END\n")


# ------------------------------ index helpers --------------------------------- #
def test_coordinating_indices_relaxed_both_vs_ddxxd_only():
    both = coordinating_indices_relaxed(_SEQ)
    assert both == _COORD_IDX
    # For a both-motif sequence, relaxed == strict both-motif set.
    assert _coordinating_indices_both(_SEQ) == both
    # DDXXD-only sequence -> relaxed keeps just the DDXXD coordinating residues.
    ddxxd_only = "AAAAA" + "DDAAD" + "AAAAAAAAAA"
    assert coordinating_indices_relaxed(ddxxd_only) == [5, 6, 9]
    assert _coordinating_indices_both(ddxxd_only) is None  # strict needs both
    # No motif -> None.
    assert coordinating_indices_relaxed("A" * 30) is None


# ------------------------------ cage geometry --------------------------------- #
def test_cage_geometry_closed_form():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "cage.pdb")
        _write_atoms(p, _cage_atoms())
        r = active_site_geometry(p)
        assert r["n_coordinating_oxygens"] == 12
        _approx(r["carboxylate_convergence_radius"], 3.0)
        # Nearest non-oxygen protein atom to the oxygen centroid (origin) is a CA at r=5.
        _approx(r["metal_point_void"], 5.0)
        assert r["n_residues"] == len(_SEQ)


def test_metal_point_is_oxygen_centroid():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "cage.pdb")
        _write_atoms(p, _cage_atoms())
        seq, residues, _ = structure_sequence_residues_atoms(p)
        mp = metal_point(seq, residues)
        assert mp is not None
        np.testing.assert_allclose(mp, [0.0, 0.0, 0.0], atol=1e-6)


def test_missing_motif_is_nan():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "poly.pdb")
        _write_poly_ala(p)
        r = active_site_geometry(p)
        assert np.isnan(r["carboxylate_convergence_radius"])
        assert np.isnan(r["metal_point_void"])
        assert r["n_coordinating_oxygens"] == 0
        assert r["best_template"] == ""


def test_ddxxd_only_still_reports_cage():
    # Relaxed set anchors on DDXXD alone -> 3 ASP * 2 oxygens = 6, not NaN.
    seq = "AAAAA" + "DDAAD" + "AAAAAAAAAA"
    # Assign the 6 ASP oxygens (indices 5,6,9) to sphere points, keyed by residue so we
    # can write atoms in RESIDUE ORDER (CA then oxygens) — otherwise the PDB residue
    # ordering, and hence the derived sequence, would be scrambled.
    slots = [(i, name) for i in (5, 6, 9) for name in ("OD1", "OD2")]
    ox_by_res = {}
    for (i, name), xyz in zip(slots, _SPHERE12):
        ox_by_res.setdefault(i, []).append((name, xyz))
    atoms = []
    for i, aa in enumerate(seq):
        atoms.append(("CA", _THREE[aa], i + 1, (100 + i, 0, 0)))
        for name, xyz in ox_by_res.get(i, []):
            atoms.append((name, "ASP", i + 1, xyz))
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "one.pdb")
        _write_atoms(p, atoms)
        r = active_site_geometry(p)
        assert r["n_coordinating_oxygens"] == 6
        assert np.isfinite(r["carboxylate_convergence_radius"])


# ------------------------------ constellation RMSD ---------------------------- #
def test_constellation_rmsd_translation_invariant():
    with tempfile.TemporaryDirectory() as d:
        _write_atoms(os.path.join(d, "ref.pdb"), _cage_atoms())
        # A rigidly shifted copy: after superposition its constellation RMSD ~ 0.
        _write_atoms(os.path.join(d, "shifted.pdb"), _cage_atoms(shift=(100.0, -30.0, 7.0)))
        templates = build_templates(d, ["ref"])
        assert "ref" in templates
        r = active_site_geometry(os.path.join(d, "shifted.pdb"), templates=templates)
        assert np.isfinite(r["catalytic_constellation_rmsd"])
        _approx(r["catalytic_constellation_rmsd"], 0.0, tol=1e-3)
        assert r["best_template"] == "ref"


# ------------------------------ dir driver ------------------------------------ #
def test_dir_id_keying_csv_and_nan_on_broken():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write_atoms(os.path.join(structs, "good.pdb"), _cage_atoms())
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("not a pdb\n")
        df = active_site_geometry_dir(structs, template_ids=["good"])

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        assert os.path.isfile(structs + "_active_site_geometry.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["carboxylate_convergence_radius"]), 3.0)
        assert int(good["n_coordinating_oxygens"]) == 12
        _approx(float(good["catalytic_constellation_rmsd"]), 0.0, tol=1e-3)

        broken = df.set_index("ID").loc["broken"]
        assert np.isnan(broken["carboxylate_convergence_radius"])
        assert int(broken["n_residues"]) == 0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
