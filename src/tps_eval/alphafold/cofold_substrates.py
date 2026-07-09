"""Canonical class-I TPS prenyl-diphosphate substrates for AlphaFold3 holo co-folding.

Single source of truth for the substrate SMILES used when ``--af3_cofold mg_<substrate>``
(force one substrate for every design) or ``--af3_cofold mg_ee`` (per-design substrate from
the EnzymeExplorer sequence-only call) co-folds the substrate alongside the trinuclear Mg2+
cluster.

SMILES are keyed by UPPERCASE code, matching the EE / ``knn.substrate_class`` substrate
vocabulary (so ``mg_ee`` and the CataPro ``--target_substrate`` parameter can map a substrate
code straight to a SMILES). ``SUBSTRATE_SMILES`` is the single source of truth for substrate
SMILES across tps_eval — the full vocabulary shared by two consumers: AF3 co-folding AND CataPro
kinetics. The clean single-molecule prenyl-diphosphates are present (DMAPP C5, GPP C10, FPP C15,
GGPP C20, GFPP C25, C35 heptaprenyl-PP). Multi-molecule / non-diphosphate EE classes (EDSQ
epoxysqualene, 2xGGPP phytoene, IDS prenyltransferase) are intentionally absent — CataPro
returns NaN for those, and ``mg_ee`` falls back to Mg-only.

Which subset AF3 actually co-folds is a SEPARATE policy: the ``COFOLDABLE`` list below
(GPP/FPP/GGPP/GFPP), NOT membership in ``SUBSTRATE_SMILES``. So adding a CataPro-only substrate
SMILES here never changes AF3 co-folding behaviour.

AF3 ligand caveat: SMILES is used (vs a PDB CCD code) for parity with the validated structures;
AF3 ligand geometry from SMILES is a hypothesis — verify the diphosphate lands at the
DDXXD/NSE cage downstream (that is exactly what the ``substrate_positioning`` tool measures).
"""
from __future__ import annotations
from typing import Dict, List

# substrate code (UPPERCASE, EE/substrate_class vocabulary) -> SMILES. Full vocabulary shared
# by AF3 co-folding and CataPro; co-foldability is governed separately by COFOLDABLE (below).
SUBSTRATE_SMILES: Dict[str, str] = {
    "DMAPP": "CC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",                                 # C5  hemi
    "GPP":  "CC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",                          # C10 mono
    "FPP":  "CC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",                  # C15 sesqui
    "GGPP": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",          # C20 di
    "GFPP": "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",  # C25 sester
    "C35":  "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O",  # C35 sesquar
}

# Co-foldable substrate codes (sorted by chain length) — the set --af3_cofold accepts as
# `mg_<code>` and the set mg_ee restricts the EE argmax to.
COFOLDABLE: List[str] = ["GPP", "FPP", "GGPP", "GFPP"]


def smiles_for(code: str) -> str:
    """SMILES for a substrate code (case-insensitive). Raises KeyError if the code has no
    known SMILES (e.g. EDSQ / 2xGGPP / IDS — CataPro treats those as NaN)."""
    return SUBSTRATE_SMILES[code.upper()]


def is_cofoldable(code: str) -> bool:
    """Whether AF3 co-folds this substrate. Governed by the explicit COFOLDABLE list, NOT by
    SUBSTRATE_SMILES membership — a substrate may have a SMILES (usable by CataPro) yet not be
    co-folded by AF3 (e.g. DMAPP, C35)."""
    return bool(code) and code.upper() in COFOLDABLE
