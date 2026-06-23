from __future__ import annotations

"""Self-contained tests for cyclization_geometry.py (numpy + biopython only — no conda
env / EnzymeExplorer / PyMOL). Run from this directory so the flat-module imports resolve:
    cd src/structure_metrics && python test_cyclization_geometry.py
or:
    cd src/structure_metrics && python -m pytest test_cyclization_geometry.py -q

Builds a tiny synthetic structure: a few ALA + one aromatic PHE near a synthetic prenyl-PP
ligand ('FPP' HETATM: 1 P + 4 O diphosphate + a 10-carbon chain). Asserts substrate
detection, the fold metrics (rgyr / fold-back / end-to-end) and the cation-pi track, plus
an apo (no-ligand) graceful not-applicable row.
"""

import os
import tempfile

import numpy as np

from cyclization_geometry import COLUMNS, cyclization_geometry, cyclization_geometry_dir

PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")


def _atom_line(serial, name, resname, chain, resseq, xyz, record="ATOM", element=None):
    x, y, z = xyz
    if element is None:
        element = name[0]
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{record:<6}{serial:>5} {atom_field}{'':1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{80.0:>6.2f}{'':10}{element:>2}\n"
    )


def _write_synthetic(path, *, with_ligand=True):
    lines = []
    serial = 1
    chain = "A"
    # A few ALA backbone residues marching along x.
    for i in range(3):
        lines.append(_atom_line(serial, "CA", "ALA", chain, i + 1, (20.0 + i * 3.8, 0, 0), element="C"))
        serial += 1
    # One PHE with a benzene ring centred near (6, 3, 0) -- close to the ligand chain so it
    # lines several substrate carbons within the cation-pi cutoff.
    ring_center = np.array([6.0, 3.0, 0.0])
    ring_offsets = np.array(
        [[0, 1.2, 0], [1.0, 0.6, 0], [1.0, -0.6, 0], [0, -1.2, 0], [-1.0, -0.6, 0], [-1.0, 0.6, 0]]
    )
    lines.append(_atom_line(serial, "CA", "PHE", chain, 4, (6.0, 6.0, 0), element="C"))
    serial += 1
    for nm, off in zip(PHE_RING, ring_offsets):
        lines.append(_atom_line(serial, nm, "PHE", chain, 4, tuple(ring_center + off), element="C"))
        serial += 1
    if with_ligand:
        # Diphosphate: one P at the origin + four O around it.
        lines.append(_atom_line(serial, "P1", "FPP", chain, 900, (0, 0, 0), record="HETATM", element="P"))
        serial += 1
        for j, o in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0)]):
            lines.append(_atom_line(serial, f"O{j+1}", "FPP", chain, 900, o, record="HETATM", element="O"))
            serial += 1
        # A 10-carbon chain along x at 1.5 A spacing (C-C bonds < 1.8 A so the BFS connects
        # them); C1 (nearest the diphosphate) leads.
        for k in range(10):
            lines.append(_atom_line(serial, f"C{k+1}", "FPP", chain, 900, (1.5 + k * 1.5, 0, 0),
                                    record="HETATM", element="C"))
            serial += 1
    lines.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(lines)


def test_substrate_fold_and_track():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "holo.pdb")
        _write_synthetic(p, with_ligand=True)
        r = cyclization_geometry(p)
        assert r["substrate_present"] is True
        assert r["n_substrate_carbons"] == 10
        assert r["substrate_rgyr"] > 0
        assert np.isfinite(r["foldback_c1_to_distal"])
        assert r["substrate_endtoend"] > 0
        assert r["n_aromatics_lining"] >= 1
        assert r["n_aromatic_carbon_contacts"] >= 1
        assert 0.0 <= r["frac_aromatic_track"] <= 1.0
        assert np.isfinite(r["mean_carbon_to_aromatic"])


def test_apo_not_applicable():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "apo.pdb")
        _write_synthetic(p, with_ligand=False)
        r = cyclization_geometry(p)
        assert r["substrate_present"] is False
        assert r["n_substrate_carbons"] == 0
        assert np.isnan(r["substrate_rgyr"])
        assert np.isnan(r["foldback_c1_to_distal"])


def test_dir_writes_csv():
    with tempfile.TemporaryDirectory() as d:
        _write_synthetic(os.path.join(d, "holo.pdb"), with_ligand=True)
        _write_synthetic(os.path.join(d, "apo.pdb"), with_ligand=False)
        out = os.path.join(d, "out.csv")
        df = cyclization_geometry_dir(d, save_path=out)
        assert os.path.exists(out)
        assert len(df) == 2
        assert list(df.columns) == COLUMNS
        assert int(df["substrate_present"].sum()) == 1


if __name__ == "__main__":
    for fn in (test_substrate_fold_and_track, test_apo_not_applicable, test_dir_writes_csv):
        fn()
        print("ok", fn.__name__)
    print("all cyclization_geometry tests passed")
