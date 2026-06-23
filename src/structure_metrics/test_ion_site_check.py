from __future__ import annotations

"""Self-contained tests for ion_site_check.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_ion_site_check.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_ion_site_check.py -q

These build tiny synthetic PDB files in a temp dir (numpy + biopython only — no
EnzymeExplorer / PyMOL / conda env): a handful of residues forming a DDXXD-like +
NSE/DTE-like motif, with the metal-coordinating side-chain oxygens placed at known
coordinates so the apo cage centroid is exactly computable, and an MG HETATM placed
either AT that centroid (well-placed) or 30 A away (mis-placed), plus an apo variant
with no HETATM. We then assert ion_in_site / well_placed / distances directly.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from ion_site_check import (
    COLUMNS,
    ion_site_check,
    ion_site_check_dir,
    read_ion_hetatms,
)

# The single-letter sequence we lay down, residue by residue. Designed so the SHARED
# motif localizer matches:
#   DDXXD-family pattern  [DE][DE]..[DE]            -> residues 0-4 "DDAAD"
#   NSE/DTE pattern       (N|D)D(L|I|V).(S|T)...E   -> residues 8-16 "NDLASGHEE"
# A few alanine spacers separate them. The coordinating residues are all
# oxygen-bearing (D/E/N/S), as metal_point requires.
#   DDXXD coordinating offsets (0,1,4) -> residues 0(ASP), 1(ASP), 4(ASP)
#   NSE   coordinating offsets (0,1,4,8) -> residues 8(ASN), 9(ASP), 12(SER), 16(GLU)
SEQUENCE = "DDAAD" + "AAA" + "NDLASGHEE" + "AAA"

# 1-letter -> 3-letter for the residue names we emit.
ONE_TO_THREE = {
    "D": "ASP",
    "E": "GLU",
    "N": "ASN",
    "S": "SER",
    "A": "ALA",
    "G": "GLY",
    "H": "HIS",
    "L": "LEU",
}

# The side-chain oxygen atom name we attach to each coordinating residue type, so the
# coordinating-oxygen helper picks it up. (active_site_geometry gathers ASP OD1/OD2,
# GLU OE1/OE2, ASN OD1, SER OG, THR OG1.)
COORD_O_ATOM = {"ASP": "OD1", "GLU": "OE1", "ASN": "OD1", "SER": "OG"}

# The 0-based residue indices whose oxygens form the cage (DDXXD + NSE coordinating).
COORD_RES_INDICES = (0, 1, 4, 8, 9, 12, 16)


def _atom_line(serial, name, resname, chain, resseq, xyz, record="ATOM", element=None):
    """One PDB ATOM/HETATM line. Biopython needs columns roughly right; we right-pad
    the atom name into cols 13-16 and use the standard fixed-width layout."""
    x, y, z = xyz
    if element is None:
        element = name[0]
    # Atom name: 4 chars. For 1-2 char names PDB convention leaves col 13 blank.
    if len(name) >= 4:
        atom_field = name[:4]
    else:
        atom_field = " " + name.ljust(3)
    return (
        f"{record:<6}{serial:>5} {atom_field}{'':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}"
        f"{'':10}{element:>2}\n"
    )


def _write_synthetic_pdb(path, *, ion_xyz=None):
    """Write a synthetic PDB with the SEQUENCE laid out along x (CA every 3.8 A), a
    side-chain coordinating oxygen on each coordinating residue clustered near the
    origin (so the cage centroid is ~origin), and optionally an MG HETATM at ion_xyz.

    Returns the cage centroid (mean of the placed coordinating oxygens)."""
    lines = []
    serial = 1
    chain = "A"

    # Coordinating oxygens: a tight cluster around the origin (within ~1 A), so their
    # centroid (the cage point) is predictable and small.
    oxygen_offsets = np.array(
        [
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
            [0.3, 0.3, 0.0],
        ]
    )
    coord_o_coords = []
    coord_o_iter = iter(oxygen_offsets)

    for i, aa in enumerate(SEQUENCE):
        resname = ONE_TO_THREE[aa]
        resseq = i + 1
        # CA marches along x, far from the origin cluster, so it never gets mistaken
        # for a coordinating oxygen or a clashing atom.
        ca_xyz = np.array([20.0 + i * 3.8, 0.0, 0.0])
        lines.append(_atom_line(serial, "CA", resname, chain, resseq, ca_xyz, element="C"))
        serial += 1
        # On coordinating residues, attach the side-chain oxygen near the origin.
        if i in COORD_RES_INDICES:
            o_xyz = next(coord_o_iter)
            coord_o_coords.append(o_xyz)
            lines.append(
                _atom_line(serial, COORD_O_ATOM[resname], resname, chain, resseq, o_xyz, element="O")
            )
            serial += 1

    centroid = np.vstack(coord_o_coords).mean(axis=0)

    if ion_xyz is not None:
        lines.append(
            _atom_line(serial, "MG", "MG", chain, 900, ion_xyz, record="HETATM", element="MG")
        )
        serial += 1

    lines.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return centroid


def test_read_ion_hetatms_selects_by_resname():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write_synthetic_pdb(p, ion_xyz=(0.0, 0.0, 0.0))
        ions, diphos = read_ion_hetatms(p)
        assert ions.shape == (1, 3)
        assert diphos.shape == (0, 3)
        # A non-matching resname filter still detects the MG by its ELEMENT (the element-based
        # fallback for AF3 'MG' / Boltz2 'LIG2'): resname filtering only ADDS resnames, it cannot
        # exclude a true Mg/Mn element atom.
        ions2, _ = read_ion_hetatms(p, ion_resnames=("ZN",))
        assert ions2.shape == (1, 3)


def test_well_placed_ion_at_centroid():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "good.pdb")
        centroid = _write_synthetic_pdb(p, ion_xyz=None)  # learn the centroid first
        # Re-write with the ion AT the cage centroid.
        _write_synthetic_pdb(p, ion_xyz=tuple(centroid))
        r = ion_site_check(p)
        assert r["metal_point_found"] is True
        assert r["n_ions_modelled"] == 1
        assert r["min_ion_to_cage_dist"] < 1e-3
        assert r["n_ions_in_site"] == 1
        assert r["ion_in_site"] is True
        # All 7 coordinating oxygens are within ~1 A of the centroid -> within the
        # 2.8 A coord cutoff -> the ion is coordinated and well-placed.
        assert r["max_coordinating_contacts"] >= 2
        assert r["n_ions_coordinated"] == 1
        assert r["well_placed"] is True
        assert r["mg_canonical_motif_coordination"] is True
        assert r["n_motif_coord_asp"] >= 1
        assert r["n_motif_coord_nse"] >= 1
        assert r["mg_to_motif_dist"] < 2.0
        assert r["n_diphosphate_atoms"] == 0
        assert np.isnan(r["diphosphate_to_cage_dist"])


def test_misplaced_ion_far_away():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "bad.pdb")
        centroid = _write_synthetic_pdb(p, ion_xyz=None)
        far = tuple(np.asarray(centroid) + np.array([30.0, 0.0, 0.0]))
        _write_synthetic_pdb(p, ion_xyz=far)
        r = ion_site_check(p)
        assert r["metal_point_found"] is True
        assert r["n_ions_modelled"] == 1
        assert r["min_ion_to_cage_dist"] > 25.0
        assert r["n_ions_in_site"] == 0
        assert r["ion_in_site"] is False
        assert r["max_coordinating_contacts"] == 0
        assert r["n_ions_coordinated"] == 0
        assert r["well_placed"] is False
        assert r["mg_canonical_motif_coordination"] is False
        assert r["n_motif_coord_asp"] == 0
        assert r["n_motif_coord_nse"] == 0


def test_apo_structure_not_applicable():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "apo.pdb")
        _write_synthetic_pdb(p, ion_xyz=None)  # no HETATM
        r = ion_site_check(p)
        # Cage is still computable (the protein has the motifs)...
        assert r["metal_point_found"] is True
        # ...but no ions -> a graceful not-applicable row.
        assert r["n_ions_modelled"] == 0
        assert np.isnan(r["min_ion_to_cage_dist"])
        assert r["n_ions_in_site"] == 0
        assert r["ion_in_site"] is False
        assert r["well_placed"] is False
        assert r["mg_canonical_motif_coordination"] is False
        assert np.isnan(r["diphosphate_to_cage_dist"])


def test_dir_driver_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        centroid = _write_synthetic_pdb(os.path.join(structs, "good.pdb"), ion_xyz=None)
        _write_synthetic_pdb(os.path.join(structs, "good.pdb"), ion_xyz=tuple(centroid))
        far = tuple(np.asarray(centroid) + np.array([30.0, 0.0, 0.0]))
        _write_synthetic_pdb(os.path.join(structs, "bad.pdb"), ion_xyz=far)
        _write_synthetic_pdb(os.path.join(structs, "apo.pdb"), ion_xyz=None)

        df = ion_site_check_dir(structs)
        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "bad", "apo"}
        assert os.path.isfile(structs + "_ion_site_check.csv")

        good = df.set_index("ID").loc["good"]
        assert bool(good["well_placed"]) is True
        assert bool(good["ion_in_site"]) is True

        bad = df.set_index("ID").loc["bad"]
        assert bool(bad["ion_in_site"]) is False
        assert float(bad["min_ion_to_cage_dist"]) > 25.0

        apo = df.set_index("ID").loc["apo"]
        assert int(apo["n_ions_modelled"]) == 0
        assert pd.isna(apo["min_ion_to_cage_dist"])
        assert bool(apo["ion_in_site"]) is False


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
