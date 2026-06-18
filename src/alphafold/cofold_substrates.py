"""Canonical class-I TPS prenyl-diphosphate substrates for AlphaFold3 holo co-folding.

Single source of truth for the substrate SMILES used when ``--af3_cofold mg_<substrate>``
(force one substrate for every design) or ``--af3_cofold mg_ee`` (per-design substrate from
the EnzymeExplorer sequence-only call) co-folds the substrate alongside the trinuclear Mg2+
cluster.

SMILES are the EnzymeExplorer substrate set (the same ligands used to produce the validated
candidate co-folds). Codes are UPPERCASE, matching the EE / ``knn.substrate_class`` substrate
vocabulary (so ``mg_ee`` can map an EE argmax substrate straight to a SMILES). Only the linear
prenyl-diphosphate substrates that have a single well-defined chain are co-foldable here;
exotic / multi-molecule EE classes (2xFPP, 2xGGPP, EDSQ epoxysqualene, CPP copalyl-PP, IDS
prenyltransferase) are intentionally excluded — a design whose EE argmax is one of those falls
back to Mg-only under ``mg_ee``.

AF3 ligand caveat: SMILES is used (vs a PDB CCD code) for parity with the validated structures;
AF3 ligand geometry from SMILES is a hypothesis — verify the diphosphate lands at the
DDXXD/NSE cage downstream (that is exactly what the ``substrate_positioning`` tool measures).
"""
from __future__ import annotations
from typing import Dict, List

# substrate code (UPPERCASE, EE/substrate_class vocabulary) -> SMILES
SUBSTRATE_SMILES: Dict[str, str] = {
    "GPP":  "CC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",                       # C10 mono
    "FPP":  "CC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",               # C15 sesqui
    "GGPP": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",       # C20 di
    "GFPP": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",  # C25 sester
}

# Co-foldable substrate codes (sorted by chain length) — the set --af3_cofold accepts as
# `mg_<code>` and the set mg_ee restricts the EE argmax to.
COFOLDABLE: List[str] = ["GPP", "FPP", "GGPP", "GFPP"]


def smiles_for(code: str) -> str:
    """SMILES for a substrate code (case-insensitive). Raises KeyError if not co-foldable."""
    return SUBSTRATE_SMILES[code.upper()]


def is_cofoldable(code: str) -> bool:
    return bool(code) and code.upper() in SUBSTRATE_SMILES
