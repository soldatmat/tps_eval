"""Unit tests for structure_conversion (mmCIF -> PDB materialization).

Regression: the vendored ProteinMPNN parses ATOM lines by fixed column offsets and
cannot read mmCIF, yet proteinmpnn_score/self_consistency handed it AF3's
``<job>_model.cif`` directly. Every structure raised, every raise became a NaN row,
and both metrics shipped an all-NaN column that read as a legitimate result.
"""
from __future__ import annotations

import os
import tempfile

from tps_eval.structure_metrics.structure_conversion import is_cif, write_pdb_copy

# Two residues of one chain plus a het (water) that must be dropped.
_CIF = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . MET A 1 1 ? 1.000 2.000 3.000 1.0 90.00 1 A 1
ATOM 2 C CA . MET A 1 1 ? 2.000 2.000 3.000 1.0 90.00 1 A 1
ATOM 3 C C . MET A 1 1 ? 3.000 2.000 3.000 1.0 90.00 1 A 1
ATOM 4 O O . MET A 1 1 ? 4.000 2.000 3.000 1.0 90.00 1 A 1
ATOM 5 N N . ALA A 1 2 ? 1.000 3.000 3.000 1.0 80.00 2 A 1
ATOM 6 C CA . ALA A 1 2 ? 2.000 3.000 3.000 1.0 80.00 2 A 1
ATOM 7 C C . ALA A 1 2 ? 3.000 3.000 3.000 1.0 80.00 2 A 1
ATOM 8 O O . ALA A 1 2 ? 4.000 3.000 3.000 1.0 80.00 2 A 1
ATOM 9 N N . GLY B 1 1 ? 8.000 2.000 3.000 1.0 70.00 1 B 1
ATOM 10 C CA . GLY B 1 1 ? 9.000 2.000 3.000 1.0 70.00 1 B 1
HETATM 11 O O . HOH C 2 1 ? 5.000 5.000 5.000 1.0 30.00 3 C 1
"""


def _write_cif(tmp: str) -> str:
    path = os.path.join(tmp, "design_model.cif")
    with open(path, "w") as fh:
        fh.write(_CIF)
    return path


def test_is_cif():
    assert is_cif("/a/b_model.cif")
    assert is_cif("/a/b.mmcif")
    assert is_cif("/a/B.CIF")
    assert not is_cif("/a/b.pdb")


def test_cif_becomes_proteinmpnn_parsable_pdb():
    """The exact failure was float(line[30:38]) on a CIF line — assert those
    columns now hold the coordinates ProteinMPNN expects."""
    tmp = tempfile.mkdtemp(prefix="structconv_")
    out = write_pdb_copy(_write_cif(tmp), os.path.join(tmp, "nested", "design.pdb"))
    assert os.path.isfile(out)
    atom_lines = [l for l in open(out) if l.startswith("ATOM")]
    assert atom_lines
    for line in atom_lines:
        x, y, z = (float(line[i:i + 8]) for i in (30, 38, 46))
        assert all(isinstance(v, float) for v in (x, y, z))
    # Waters (HETATM) are dropped; both polymer chains survive by default.
    assert not any(l.startswith("HETATM") for l in open(out))
    chains = {l[21] for l in atom_lines}
    assert chains == {"A", "B"}


def test_chain_id_restricts_output():
    tmp = tempfile.mkdtemp(prefix="structconv_chain_")
    out = write_pdb_copy(_write_cif(tmp), os.path.join(tmp, "a.pdb"), chain_id="A")
    chains = {l[21] for l in open(out) if l.startswith("ATOM")}
    assert chains == {"A"}


def test_pdb_input_roundtrips():
    tmp = tempfile.mkdtemp(prefix="structconv_pdb_")
    first = write_pdb_copy(_write_cif(tmp), os.path.join(tmp, "one.pdb"))
    second = write_pdb_copy(first, os.path.join(tmp, "two.pdb"))
    a = [l[:54] for l in open(first) if l.startswith("ATOM")]
    b = [l[:54] for l in open(second) if l.startswith("ATOM")]
    assert a == b
