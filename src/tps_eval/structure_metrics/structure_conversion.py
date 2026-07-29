"""Materialize a PDB copy of a structure — for tools that cannot read mmCIF.

The vendored ProteinMPNN parses coordinates by FIXED COLUMN OFFSETS out of lines
starting with ``ATOM`` (``protein_mpnn_utils.parse_PDB_biounits``), so handing it an
mmCIF does not "just work": it slices the wrong characters and dies with
``ValueError: could not convert string to float: ' 5   ? -'``. Our AF3 runs produce
``af_output/<job>/<job>_model.cif``, i.e. exactly the layout that hits this — which is
why `proteinmpnn_score` and `self_consistency` wrote all-NaN columns on AF3 inputs
(each per-structure failure was caught and turned into a NaN row, so the tools
"succeeded").

Converting through Biopython also normalizes two other things ProteinMPNN is picky
about: only the first model is kept, and only standard polymer residues are written.
"""

from __future__ import annotations

import os
import warnings
from typing import Optional

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)


def is_cif(path: str) -> bool:
    return path.lower().endswith((".cif", ".mmcif"))


def parser_for(path: str):
    """Biopython parser matching the file extension (mirrors plddt.py)."""
    return _CIF_PARSER if is_cif(path) else _PDB_PARSER


class _PolymerSelect(Select):
    """Standard polymer residues, optionally restricted to one chain."""

    def __init__(self, chain_id: Optional[str] = None):
        self.chain_id = chain_id

    def accept_chain(self, chain):  # noqa: N802 (Biopython API)
        return 1 if (self.chain_id is None or chain.id == self.chain_id) else 0

    def accept_residue(self, residue):  # noqa: N802
        return 1 if residue.id[0] == " " else 0


def write_pdb_copy(src_path: str, out_path: str, *, chain_id: Optional[str] = None) -> str:
    """Write `src_path` (.pdb or .cif) as a PDB at `out_path`; return `out_path`.

    Keeps the first model and standard polymer residues only. `chain_id` restricts
    the output to a single chain (used for monomer designability on multimers).
    """
    parser = parser_for(src_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", src_path)
    for model in list(structure)[1:]:
        structure.detach_child(model.id)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, _PolymerSelect(chain_id))
    return out_path
