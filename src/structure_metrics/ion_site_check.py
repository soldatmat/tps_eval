"""Ion-placement (catalytic-site) check for AF3 holo-folded class-I TPS designs.

Class-I terpene synthases bind a trinuclear Mg2+ (or Mn2+) cluster chelated by the
two metal-binding motifs (the DDXXD "aspartate-rich" + the NSE/DTE triad). The other
active-site structure tools are *apo-robust*: they anchor on the
``active_site_geometry.metal_point`` — the centroid of the coordinating side-chain
oxygens of those motifs — which is a HYPOTHETICAL metal locus derived from the
protein alone, NOT an actual modelled ion. AlphaFold3 can additionally co-fold a
design WITH the ions (the ``--af3_cofold mg|mg_ppi`` feature places CCD ``MG`` x3 and
optionally a ``POP`` pyrophosphate), but AF3 free-ion placement is just a hypothesis:
nothing checks whether the predicted ions actually land in the catalytic cage.

This tool is that geometric check. Unlike every other structure tool (which SKIP the
HETATM records), it READS the ion HETATMs and compares them to the expected apo cage
point:

For each structure in ``structs_dir`` we

1. parse it with Biopython (REUSING the af_output/flat auto-detection + ID-stem
   convention via ``plddt._collect_structures`` and the sequence/residue derivation
   from ``active_site_geometry.structure_sequence_residues_atoms``),
2. derive the expected apo **cage point** via the canonical, single-source-of-truth
   ``active_site_geometry.metal_point`` (DDXXD required, NSE/DTE added when matched);
   it is None when DDXXD is absent,
3. read the **ion HETATMs** (a helper that, unlike the apo parser, KEEPS HETATMs and
   selects by residue/atom name — default ion resnames ``{MG, MN}``; a separate
   diphosphate set ``{POP, PPV, PPK}`` for the ``mg_ppi`` case), and
4. measure how well the modelled ions sit in the apo cage (distance to the centroid,
   how many fall inside the site sphere, and how many cage carboxylate oxygens each
   ion coordinates at Mg-O bonding distance).

Columns (keyed by ``ID``):

* ``metal_point_found``        — bool; True when the apo cage centroid is computable
  (DDXXD present + its coordinating oxygens found). A real RED FLAG for a class-I
  design when False — recorded, never silently dropped.
* ``n_ions_modelled``          — count of ion HETATMs found (0 for apo / ESMFold).
* ``min_ion_to_cage_dist``     — nearest ion -> cage centroid distance (A); NaN when
  there are no ions OR no metal point.
* ``n_ions_in_site``           — ions within ``site_radius`` (default 5.0 A) of the cage.
* ``ion_in_site`` (bool)       — >=1 ion within ``site_radius``.
* ``max_coordinating_contacts`` — max over ions of the number of cage carboxylate/
  hydroxyl oxygens within ``coord_cutoff`` (default 2.8 A; real Mg-O ~2.0-2.5 A).
* ``n_ions_coordinated``       — ions with >=``min_coord_contacts`` coordinating-O contacts.
* ``well_placed`` (bool)       — >=1 ion coordinated by >=``min_coord_contacts`` cage
  oxygens (the strict validation that an ion sits in the carboxylate cage).
* ``n_diphosphate_atoms``      — count of diphosphate HETATM atoms (the ``mg_ppi`` case).
* ``diphosphate_to_cage_dist`` — diphosphate-atom centroid -> cage centroid distance
  (A); NaN when no diphosphate atoms OR no metal point.
* ``n_residues``               — modelled residue count, for context.

SEMANTICS / NOT-APPLICABLE: apo or ESMFold structures carry no ions, so
``n_ions_modelled=0``, the distance columns are NaN and the bool columns False — a
graceful "not applicable" row. The tool still produces a complete CSV and exits 0 in
that case (the output IS created). Only a genuine fatal error (bad dir, unreadable
input) raises / exits non-zero. This tool therefore only carries signal for AF3 holo
folds (``--af3_cofold mg|mg_ppi``).

Citation: Christianson, D. W. "Structural and Chemical Biology of Terpenoid
Cyclases." Chem. Rev. 2017, 117, 11570-11648 (the trinuclear Mg2+ cluster geometry).
"""

from __future__ import annotations

import argparse
import os
import warnings
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Reuse the shared structure loader + the carboxylate-cage machinery (single source
# of truth for the metal point) from the neighbouring structure tools. Do NOT
# re-encode the motifs or the structure parser beyond the ion-reader addition below.
from active_site_geometry import (
    COORDINATING_OXYGEN_ATOMS,
    _coordinating_oxygens,
    coordinating_indices_relaxed,
    metal_point as _cage_metal_point,
    structure_sequence_residues_atoms,
    _parser_for,
)
from plddt import _collect_structures

# --------------------------------------------------------------------------- #
# Defaults (CLI-overridable)
# --------------------------------------------------------------------------- #
# Ion HETATM residue names to read. CCD MG is what AF3 places for the trinuclear
# cluster; MN covers the Mn2+-using class-I TPS. Configurable.
DEFAULT_ION_RESNAMES = ("MG", "MN")

# Metal-ion ELEMENTS. Ion detection is element-based (not just resname) so it works across
# structure sources: AF3 names ions by resname 'MG'; Boltz2 names them 'LIG2' (element still
# MG). A monatomic residue whose atoms are all these elements is an ion regardless of resname.
ION_ELEMENTS = {"MG", "MN"}

# Diphosphate / pyrophosphate HETATM residue names (the mg_ppi case). CCD POP is the
# pyrophosphate AF3 places; PPV/PPK are alternative diphosphate ligand codes.
DEFAULT_DIPHOSPHATE_RESNAMES = ("POP", "PPV", "PPK")

# An ion within this distance (A) of the apo cage centroid counts as "in the site".
# ~5 A is generous: the cage centroid is the oxygen centroid, and the three Mg sit a
# couple A off it; a well-placed cluster lands inside this sphere.
DEFAULT_SITE_RADIUS = 5.0

# Mg-O coordination distance cutoff (A). Real octahedral Mg-O bonds are ~2.0-2.5 A;
# 2.8 A is a slightly relaxed cutoff (AF3 ion placement is not bond-length exact).
DEFAULT_COORD_CUTOFF = 2.8

# An ion counts as "coordinated" / "well placed" when it contacts at least this many
# cage carboxylate/hydroxyl oxygens within coord_cutoff.
DEFAULT_MIN_COORD_CONTACTS = 2

COLUMNS = [
    "ID",
    "metal_point_found",
    "n_ions_modelled",
    "min_ion_to_cage_dist",
    "n_ions_in_site",
    "ion_in_site",
    "max_coordinating_contacts",
    "n_ions_coordinated",
    "well_placed",
    "n_diphosphate_atoms",
    "diphosphate_to_cage_dist",
    "n_residues",
]


def read_ion_hetatms(
    structure_path: str,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    diphosphate_resnames: Tuple[str, ...] = DEFAULT_DIPHOSPHATE_RESNAMES,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read the ION HETATM atoms from a structure.

    Unlike ``structure_sequence_residues_atoms`` (which deliberately SKIPS HETATMs),
    this KEEPS them and selects by residue name: returns
    ``(ion_coords, diphosphate_coords)`` as ``(n, 3)`` arrays.

    * ``ion_coords`` — one row per ion HETATM atom whose residue name is in
      ``ion_resnames`` (monatomic ions: one atom == one ion).
    * ``diphosphate_coords`` — every atom of any HETATM residue whose name is in
      ``diphosphate_resnames`` (a diphosphate is polyatomic).

    Residue-name matching is case-insensitive and whitespace-stripped (PDB pads atom/
    residue names). Only the first model is read (predicted structures write one)."""
    ion_set = {r.strip().upper() for r in ion_resnames}
    diphos_set = {r.strip().upper() for r in diphosphate_resnames}

    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))  # first model only

    ions: List[np.ndarray] = []
    diphos: List[np.ndarray] = []
    for chain in model:
        for residue in chain:
            # HETATM only (hetflag != " "); protein residues are handled by the apo
            # parser elsewhere.
            if residue.id[0] == " ":
                continue
            resname = residue.get_resname().strip().upper()
            atom_els = [(a.element or "").strip().upper() for a in residue]
            if (resname in ion_set) or (atom_els and all(el in ION_ELEMENTS for el in atom_els)):
                for atom in residue:                      # element-based: AF3 'MG' & Boltz2 'LIG2'
                    ions.append(np.asarray(atom.get_coord(), dtype=float))
            elif resname in diphos_set:
                for atom in residue:
                    diphos.append(np.asarray(atom.get_coord(), dtype=float))
    ion_arr = np.vstack(ions) if ions else np.empty((0, 3))
    diphos_arr = np.vstack(diphos) if diphos else np.empty((0, 3))
    return ion_arr, diphos_arr


def ion_site_check(
    structure_path: str,
    *,
    site_radius: float = DEFAULT_SITE_RADIUS,
    coord_cutoff: float = DEFAULT_COORD_CUTOFF,
    min_coord_contacts: int = DEFAULT_MIN_COORD_CONTACTS,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    diphosphate_resnames: Tuple[str, ...] = DEFAULT_DIPHOSPHATE_RESNAMES,
) -> Dict[str, object]:
    """Ion-placement metrics for one structure. Apo / ESMFold (no ions) -> a graceful
    not-applicable row (n_ions_modelled=0, distances NaN, bools False)."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, object] = {
        "metal_point_found": False,
        "n_ions_modelled": 0,
        "min_ion_to_cage_dist": np.nan,
        "n_ions_in_site": 0,
        "ion_in_site": False,
        "max_coordinating_contacts": 0,
        "n_ions_coordinated": 0,
        "well_placed": False,
        "n_diphosphate_atoms": 0,
        "diphosphate_to_cage_dist": np.nan,
        "n_residues": len(sequence),
    }

    ion_coords, diphos_coords = read_ion_hetatms(
        structure_path, ion_resnames=ion_resnames, diphosphate_resnames=diphosphate_resnames
    )
    result["n_ions_modelled"] = int(len(ion_coords))
    result["n_diphosphate_atoms"] = int(len(diphos_coords))

    # Expected apo cage point (canonical relaxed coordinating-oxygen centroid).
    cage = _cage_metal_point(sequence, residues)
    if cage is None:
        # No anchor: distances to the cage are undefined (stay NaN). Ion counts are
        # still recorded above (a holo structure with no DDXXD is itself a red flag).
        return result
    result["metal_point_found"] = True

    # Gather the cage coordinating oxygens once (for the Mg-O contact count). This is
    # the SAME relaxed set the centroid was built from.
    idx = coordinating_indices_relaxed(sequence)
    cage_oxygens = _coordinating_oxygens(idx, residues) if idx is not None else np.empty((0, 3))

    if len(ion_coords):
        d_to_cage = np.sqrt(((ion_coords - cage) ** 2).sum(axis=1))
        result["min_ion_to_cage_dist"] = float(d_to_cage.min())
        n_in_site = int((d_to_cage <= site_radius).sum())
        result["n_ions_in_site"] = n_in_site
        result["ion_in_site"] = bool(n_in_site > 0)

        # Per-ion coordinating-oxygen contacts (Mg-O within coord_cutoff).
        if len(cage_oxygens):
            cutoff2 = coord_cutoff * coord_cutoff
            contacts = np.array(
                [int((((cage_oxygens - ion) ** 2).sum(axis=1) <= cutoff2).sum()) for ion in ion_coords]
            )
            result["max_coordinating_contacts"] = int(contacts.max())
            n_coord = int((contacts >= min_coord_contacts).sum())
            result["n_ions_coordinated"] = n_coord
            result["well_placed"] = bool(n_coord > 0)

    # Diphosphate centroid -> cage distance (mg_ppi case).
    if len(diphos_coords):
        diphos_centroid = diphos_coords.mean(axis=0)
        result["diphosphate_to_cage_dist"] = float(np.sqrt(((diphos_centroid - cage) ** 2).sum()))

    return result


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_ion_site_check.csv")


def ion_site_check_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    site_radius: float = DEFAULT_SITE_RADIUS,
    coord_cutoff: float = DEFAULT_COORD_CUTOFF,
    min_coord_contacts: int = DEFAULT_MIN_COORD_CONTACTS,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    diphosphate_resnames: Tuple[str, ...] = DEFAULT_DIPHOSPHATE_RESNAMES,
) -> pd.DataFrame:
    """Ion-placement check for every structure in ``structs_dir``; CSV keyed by ID.

    Mirrors the other structure-branch tools: auto-detects an AF3 ``af_output`` dir vs
    a flat dir of .pdb/.cif (via ``plddt._collect_structures``), writes
    ``<structs_dir>_ion_site_check.csv`` by default, one row per structure. Apo /
    ESMFold structures (no ions) get a graceful not-applicable row — the CSV is still
    written in full."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")
    print(
        f"Site radius {site_radius} A; coordination cutoff {coord_cutoff} A "
        f"(>= {min_coord_contacts} contacts); ions {sorted({r.upper() for r in ion_resnames})}; "
        f"diphosphate {sorted({r.upper() for r in diphosphate_resnames})}."
    )

    rows: List[Dict[str, object]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = ion_site_check(
                path,
                site_radius=site_radius,
                coord_cutoff=coord_cutoff,
                min_coord_contacts=min_coord_contacts,
                ion_resnames=ion_resnames,
                diphosphate_resnames=diphosphate_resnames,
            )
        except Exception as exc:  # malformed/unparsable -> not-applicable row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {
                "metal_point_found": False,
                "n_ions_modelled": 0,
                "min_ion_to_cage_dist": np.nan,
                "n_ions_in_site": 0,
                "ion_in_site": False,
                "max_coordinating_contacts": 0,
                "n_ions_coordinated": 0,
                "well_placed": False,
                "n_diphosphate_atoms": 0,
                "diphosphate_to_cage_dist": np.nan,
                "n_residues": 0,
            }
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    n_with_ions = int((df["n_ions_modelled"] > 0).sum())
    n_well = int(df["well_placed"].sum())
    print(
        f"Wrote {len(df)} rows to {save_path} "
        f"({n_with_ions}/{len(df)} carry modelled ions, {n_well} well-placed"
        + (f", {n_failed} unparsable" if n_failed else "")
        + ")."
    )
    if n_with_ions == 0:
        print(
            "  [note] no structures carried modelled ions -> all rows not-applicable "
            "(this tool only carries signal for AF3 holo folds: --af3_cofold mg|mg_ppi)."
        )
    return df
