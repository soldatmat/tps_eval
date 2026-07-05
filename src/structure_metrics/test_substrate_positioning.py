from __future__ import annotations

"""Self-contained tests for substrate_positioning.py (AF3 holo catalytic-site check).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_substrate_positioning.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_substrate_positioning.py -q

No external tool runs. We build tiny synthetic PDBs (numpy + biopython only) laying down a
DDXXD + NSE/DTE-motif protein whose 7 metal-coordinating side-chain oxygens sit at KNOWN
coordinates (so the carboxylate-cage centroid is exactly (0,0,0)), plus a prenyl-PP HETATM
ligand ('GPP': a P+4O diphosphate cluster and 5 chain carbons) and an MG ion at known
positions. All reported distances (diphosphate->cage, min-diphosphate->cage-oxygen,
diphosphate->ion, reactive-carbon->cage/ion) are then checked against an INDEPENDENT numpy
recomputation from the same placed coordinates, so the closed-form geometry is exact.

Covered: ligand auto-detection by composition (>=1 P and >= min_carbons C), the min_carbons
threshold, forced --substrate_resname, ion element detection, the not-applicable apo row
(no ligand -> substrate_present False, geometry NaN), and the dir driver (ID keying, CSV
naming, column order, sort, graceful failure row).

The single DDXXD/NSE motif layout mirrors test_ion_site_check.py's, which is the vetted way
to make active_site_geometry.metal_point / coordinating_indices_relaxed resolve on synthetic
input.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from substrate_positioning import (
    COLUMNS,
    read_substrate_ligand,
    substrate_positioning,
    substrate_positioning_dir,
    _default_save_path,
)

# Same motif-bearing sequence as the ion_site_check test: DDXXD family at 0-4,
# NSE/DTE at 8-16. Coordinating residues (oxygen-bearing) at 0,1,4,8,9,12,16.
SEQUENCE = "DDAAD" + "AAA" + "NDLASGHEE" + "AAA"
ONE_TO_THREE = {"D": "ASP", "E": "GLU", "N": "ASN", "S": "SER",
                "A": "ALA", "G": "GLY", "H": "HIS", "L": "LEU"}
COORD_O_ATOM = {"ASP": "OD1", "GLU": "OE1", "ASN": "OD1", "SER": "OG"}
COORD_RES_INDICES = (0, 1, 4, 8, 9, 12, 16)

# 7 coordinating-oxygen positions whose mean is exactly the origin -> cage centroid (0,0,0).
CAGE_OXYGENS = np.array([
    [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0], [0.0, -0.5, 0.0],
    [0.0, 0.0, 0.5], [0.0, 0.0, -0.5],
    [0.0, 0.0, 0.0],
])


def _atom_line(serial, name, resname, chain, resseq, xyz, record="ATOM", element=None):
    x, y, z = xyz
    if element is None:
        element = name[0]
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{record:<6}{serial:>5} {atom_field}{'':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}"
        f"{'':10}{element:>2}\n"
    )


def _write_pdb(path, *, diphos=None, carbons=None, ion=None, lig_resname="GPP",
               lig_bfactor=77.0):
    """Write the motif protein (CA along x, coordinating O at CAGE_OXYGENS) plus an
    optional prenyl-PP ligand (P + O = `diphos`; C = `carbons`) and an MG `ion`.

    Returns dict with the placed arrays for independent expected-value recomputation.
    """
    lines = []
    serial = 1
    chain = "A"
    ox_iter = iter(CAGE_OXYGENS)
    for i, aa in enumerate(SEQUENCE):
        resname = ONE_TO_THREE[aa]
        resseq = i + 1
        ca_xyz = np.array([50.0 + i * 3.8, 0.0, 0.0])   # far from origin cage
        lines.append(_atom_line(serial, "CA", resname, chain, resseq, ca_xyz, element="C"))
        serial += 1
        if i in COORD_RES_INDICES:
            o_xyz = next(ox_iter)
            lines.append(_atom_line(serial, COORD_O_ATOM[resname], resname, chain, resseq,
                                    o_xyz, element="O"))
            serial += 1

    # Prenyl-PP ligand: diphosphate (1 P + 4 O), then chain carbons. Shared resseq 900.
    if diphos is not None:
        p, o1, o2, o3, o4 = diphos
        for name, el, xyz in [("P1", "P", p), ("O1", "O", o1), ("O2", "O", o2),
                              ("O3", "O", o3), ("O4", "O", o4)]:
            line = _atom_line(serial, name, lig_resname, chain, 900, xyz,
                              record="HETATM", element=el)
            # override the fixed 80.00 b-factor with lig_bfactor
            line = line[:60] + f"{lig_bfactor:>6.2f}" + line[66:]
            lines.append(line)
            serial += 1
    if carbons is not None:
        for k, cxyz in enumerate(carbons, start=1):
            line = _atom_line(serial, f"C{k}", lig_resname, chain, 900, cxyz,
                              record="HETATM", element="C")
            line = line[:60] + f"{lig_bfactor:>6.2f}" + line[66:]
            lines.append(line)
            serial += 1
    if ion is not None:
        lines.append(_atom_line(serial, "MG", "MG", chain, 901, ion,
                                record="HETATM", element="MG"))
        serial += 1

    lines.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return {
        "cage": CAGE_OXYGENS.mean(axis=0),
        "cage_oxygens": CAGE_OXYGENS,
        "diphos": None if diphos is None else np.asarray(diphos, float),
        "carbons": None if carbons is None else np.asarray(carbons, float),
        "ion": None if ion is None else np.asarray([ion], float),
    }


# A well-placed substrate: diphosphate cluster centered at (2,0,0); 5 carbons marching
# out along +x (nearest to the diphosphate is C1 at (3,0,0)); Mg at (1,0,0).
DIPHOS = [(2.0, 0, 0), (2.5, 0, 0), (1.5, 0, 0), (2.0, 0.5, 0), (2.0, -0.5, 0)]
CARBONS = [(3.0, 0, 0), (4.0, 0, 0), (5.0, 0, 0), (6.0, 0, 0), (7.0, 0, 0)]
ION = (1.0, 0, 0)


def _approx(a, b, tol=1e-4):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_read_substrate_ligand_autodetect():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write_pdb(p, diphos=DIPHOS, carbons=CARBONS, ion=ION)
        resname, all_c, all_b, diphos, carbons, ions = read_substrate_ligand(p)
        assert resname == "GPP"
        assert all_c.shape == (10, 3)      # 5 diphosphate + 5 carbon
        assert diphos.shape == (5, 3)
        assert carbons.shape == (5, 3)
        assert ions.shape == (1, 3)
        _approx(float(all_b.mean()), 77.0, tol=1e-2)


def test_read_substrate_ligand_min_carbons_threshold():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write_pdb(p, diphos=DIPHOS, carbons=CARBONS, ion=ION)
        # Raise the carbon requirement above the ligand's 5 C -> not auto-detected.
        resname, all_c, _b, _dp, _c, ions = read_substrate_ligand(p, min_carbons=6)
        assert resname == ""
        assert all_c.shape == (0, 3)
        assert ions.shape == (1, 3)        # ion still found (element-based)
        # ...but forcing the resname detects it regardless of the carbon count.
        resname2, all_c2, _b2, _dp2, _c2, _i2 = read_substrate_ligand(
            p, min_carbons=6, substrate_resname="GPP")
        assert resname2 == "GPP"
        assert all_c2.shape == (10, 3)


def test_substrate_positioning_geometry_closed_form():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        placed = _write_pdb(p, diphos=DIPHOS, carbons=CARBONS, ion=ION)
        r = substrate_positioning(p)

        cage = placed["cage"]                      # (0,0,0)
        diphos = placed["diphos"]
        carbons = placed["carbons"]
        ion = placed["ion"]
        cage_ox = placed["cage_oxygens"]

        assert r["substrate_present"] is True
        assert r["substrate_resname"] == "GPP"
        assert r["metal_point_found"] is True
        assert int(r["n_substrate_atoms"]) == 10
        _approx(float(r["substrate_plddt"]), 77.0, tol=1e-2)
        assert int(r["n_residues"]) == len(SEQUENCE)

        # Independent numpy recomputation of every geometric column.
        exp_dp_cage = float(np.linalg.norm(diphos.mean(0) - cage))
        _approx(r["diphosphate_to_cage_dist"], exp_dp_cage)
        _approx(r["diphosphate_to_cage_dist"], 2.0)

        exp_min_o = float(np.sqrt(((diphos[:, None, :] - cage_ox[None, :, :]) ** 2)
                                  .sum(2)).min())
        _approx(r["min_diphosphate_to_cage_oxygen"], exp_min_o)
        _approx(r["min_diphosphate_to_cage_oxygen"], 1.0)
        assert r["substrate_in_site"] is True      # 1.0 <= coord_cutoff (4.0)

        exp_dp_ion = float(np.sqrt(((diphos[:, None, :] - ion[None, :, :]) ** 2)
                                   .sum(2)).min())
        _approx(r["diphosphate_to_nearest_ion"], exp_dp_ion)
        _approx(r["diphosphate_to_nearest_ion"], 0.5)

        _approx(r["diphosphate_to_ion_centroid"],
                float(np.linalg.norm(diphos.mean(0) - ion.mean(0))))
        _approx(r["diphosphate_to_ion_centroid"], 1.0)

        # Reactive carbon = the carbon nearest the diphosphate (C1 at (3,0,0)).
        d_c = np.sqrt(((carbons[:, None, :] - diphos[None, :, :]) ** 2).sum(2)).min(1)
        reactive_c = carbons[int(np.argmin(d_c))]
        _approx(r["reactive_carbon_to_cage_dist"], float(np.linalg.norm(reactive_c - cage)))
        _approx(r["reactive_carbon_to_cage_dist"], 3.0)
        _approx(r["reactive_carbon_to_nearest_ion"],
                float(np.sqrt(((reactive_c[None, :] - ion) ** 2).sum(1)).min()))
        _approx(r["reactive_carbon_to_nearest_ion"], 2.0)
        _approx(r["reactive_carbon_to_ion_centroid"],
                float(np.linalg.norm(reactive_c - ion.mean(0))))


def test_substrate_far_from_cage_not_in_site():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "far.pdb")
        far_diphos = [(x + 30.0, y, z) for (x, y, z) in DIPHOS]
        far_carbons = [(x + 30.0, y, z) for (x, y, z) in CARBONS]
        _write_pdb(p, diphos=far_diphos, carbons=far_carbons, ion=ION)
        r = substrate_positioning(p)
        assert r["substrate_present"] is True
        assert r["metal_point_found"] is True
        assert float(r["diphosphate_to_cage_dist"]) > 25.0
        assert float(r["min_diphosphate_to_cage_oxygen"]) > 25.0
        assert r["substrate_in_site"] is False


def test_apo_not_applicable():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "apo.pdb")
        _write_pdb(p)   # no ligand, no ion
        r = substrate_positioning(p)
        assert r["substrate_present"] is False
        assert r["substrate_resname"] == ""
        assert int(r["n_substrate_atoms"]) == 0
        # No substrate -> geometry NaN and metal_point_found stays False (early return).
        assert r["metal_point_found"] is False
        for k in ("diphosphate_to_cage_dist", "min_diphosphate_to_cage_oxygen",
                  "diphosphate_to_nearest_ion", "reactive_carbon_to_cage_dist",
                  "substrate_plddt"):
            assert np.isnan(r[k]), k
        assert r["substrate_in_site"] is False
        assert int(r["n_residues"]) == len(SEQUENCE)


def test_default_save_path():
    assert _default_save_path("/a/b/structs") == \
        os.path.join("/a/b", "structs_substrate_positioning.csv")


def test_dir_driver_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write_pdb(os.path.join(structs, "holo.pdb"),
                   diphos=DIPHOS, carbons=CARBONS, ion=ION)
        _write_pdb(os.path.join(structs, "apo.pdb"))
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("this is not a structure\n")

        df = substrate_positioning_dir(structs)
        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"holo", "apo", "broken"}
        assert list(df["ID"]) == ["apo", "broken", "holo"]   # sorted
        assert os.path.isfile(structs + "_substrate_positioning.csv")

        holo = df.set_index("ID").loc["holo"]
        assert bool(holo["substrate_present"]) is True
        assert bool(holo["substrate_in_site"]) is True
        _approx(float(holo["diphosphate_to_cage_dist"]), 2.0)

        apo = df.set_index("ID").loc["apo"]
        assert bool(apo["substrate_present"]) is False
        assert pd.isna(apo["diphosphate_to_cage_dist"])

        broken = df.set_index("ID").loc["broken"]     # graceful failure row
        assert bool(broken["substrate_present"]) is False
        assert int(broken["n_residues"]) == 0
        assert pd.isna(broken["diphosphate_to_cage_dist"])


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
