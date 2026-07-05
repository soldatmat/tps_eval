from __future__ import annotations

"""Self-contained tests for extract_pdb_files.py (AF3 mmCIF -> PDB conversion).

Run from this directory:
    cd src/alphafold && python test_extract_pdb_files.py
or under pytest:
    cd src/alphafold && python -m pytest test_extract_pdb_files.py -q

Uses only Biopython + tempfile (no AF3). We build a tiny structure in memory (one
protein residue, one long-resName ligand 'LIG_B', one water), write it out as mmCIF
via MMCIFIO, then run cif_to_pdb_sanitized and reparse the PDB to assert: the >3-char
resName is truncated to a valid 3-char code, the ligand HETATM is kept, water is
dropped, and the B-factor (pLDDT) survives the round-trip. Helpers _ResidueSelect and
_sanitize_resnames are also unit-tested directly.
"""

import os
import sys
import tempfile
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from alphafold.extract_pdb_files import (
    _ResidueSelect,
    _sanitize_resnames,
    cif_to_pdb_sanitized,
)


def _atom(name, coord, bfactor, element):
    return Atom(name, coord, bfactor, 1.0, " ", name, 1, element=element)


def _build_structure():
    """protein ALA (chain A, res 1) + ligand 'LIG_B' HETATM (res 2) + water (res 3)."""
    s = Structure("s")
    m = Model(0)
    c = Chain("A")
    prot = Residue((" ", 1, " "), "ALA", "")
    prot.add(_atom("CA", (0.0, 0.0, 0.0), 88.5, "C"))
    lig = Residue(("H_LIG_B", 2, " "), "LIG_B", "")
    lig.add(_atom("C1", (5.0, 0.0, 0.0), 50.0, "C"))
    wat = Residue(("W", 3, " "), "HOH", "")
    wat.add(_atom("O", (9.0, 9.0, 9.0), 10.0, "O"))
    for r in (prot, lig, wat):
        c.add(r)
    m.add(c)
    s.add(m)
    return s


def test_sanitize_resnames_truncates_long():
    s = _build_structure()
    _sanitize_resnames(s)
    names = {r.resname for r in s.get_residues()}
    assert "LIG" in names, names           # LIG_B -> LIG
    assert "LIG_B" not in names
    assert "ALA" in names and "HOH" in names  # <=3-char names untouched


def test_residue_select_keeps_hetero_drops_water():
    sel = _ResidueSelect(keep_hetero=True)
    s = _build_structure()
    res = {r.resname: r for r in s.get_residues()}
    assert sel.accept_residue(res["ALA"]) is True
    assert sel.accept_residue(res["HOH"]) is False
    assert sel.accept_residue(res["LIG_B"]) is True
    # keep_hetero=False drops the ligand too, keeps the protein.
    sel_no = _ResidueSelect(keep_hetero=False)
    assert sel_no.accept_residue(res["ALA"]) is True
    assert sel_no.accept_residue(res["LIG_B"]) is False


def test_cif_to_pdb_roundtrip_keeps_ligand_bfactor_drops_water():
    with tempfile.TemporaryDirectory() as d:
        cif = os.path.join(d, "job_model.cif")
        pdb = os.path.join(d, "job.pdb")
        io = MMCIFIO()
        io.set_structure(_build_structure())
        io.save(cif)

        cif_to_pdb_sanitized(cif, pdb, keep_hetero=True)

        out = PDBParser(QUIET=True).get_structure("o", pdb)
        resnames = [r.resname.strip() for r in out.get_residues()]
        assert "ALA" in resnames
        assert "LIG" in resnames         # truncated, parseable
        assert "LIG_B" not in resnames
        assert "HOH" not in resnames     # water dropped

        # B-factor (pLDDT) preserved on the protein CA.
        ca = next(a for a in out.get_atoms() if a.get_name() == "CA")
        assert abs(ca.get_bfactor() - 88.5) < 1e-2, ca.get_bfactor()


def test_cif_to_pdb_no_hetero_drops_ligand():
    with tempfile.TemporaryDirectory() as d:
        cif = os.path.join(d, "job_model.cif")
        pdb = os.path.join(d, "job.pdb")
        io = MMCIFIO()
        io.set_structure(_build_structure())
        io.save(cif)
        cif_to_pdb_sanitized(cif, pdb, keep_hetero=False)
        out = PDBParser(QUIET=True).get_structure("o", pdb)
        resnames = [r.resname.strip() for r in out.get_residues()]
        assert resnames == ["ALA"], resnames


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
