from __future__ import annotations

"""Cyclization-relevant holo geometry for AF3-cofolded class-I TPS designs.

Where ``ion_site_check`` / ``substrate_positioning`` confirm that the metal cluster is
assembled and the prenyl-PP substrate is BOUND and poised to ionize, this tool reports
two cheap, *necessary-not-sufficient* geometric signals for whether the bound substrate
is organized to CYCLIZE:

  (1) SUBSTRATE FOLD - is the prenyl chain curled so a distal carbon can reach C1 (the
      first ring-closure geometry), vs splayed into an extended/unproductive pose? ->
      substrate radius-of-gyration, C1->distal-carbon fold-back distance, end-to-end span.
  (2) CATION-pi TRACK - are aromatic side chains lined along the substrate carbons to
      stabilize the migrating carbocation cascade? -> number/fraction of substrate carbons
      within cation-pi range of an aromatic ring centroid, and the count of lining aromatics.

It READS the prenyl-PP ligand HETATMs (reusing ``substrate_positioning.read_substrate_ligand``,
which is validated across AF3 'MG'/FPP and Boltz2 'LIG2' naming) and the protein aromatic
residues. Fold-agnostic (.pdb/.cif). Apo / ESMFold / Mg-only / mg_ppi structures carry no
prenyl-PP substrate -> a graceful not-applicable row (substrate_present=False, geometry NaN).

These metrics are REFERENCE-INDEPENDENT (they never use the apo metal_point): the fold metrics
are intrinsic to the ligand, the cation-pi track is substrate-carbon -> protein-aromatic. They
are NECESSARY-NOT-SUFFICIENT for cyclization - a folded, aromatic-lined substrate is consistent
with a competent cyclase but does NOT establish the product (the carbocation cascade, hydride/
methyl shifts and quench are not modelled).

Columns (keyed by ``ID``):
* ``substrate_present``            - bool; a prenyl-PP substrate ligand was found.
* ``n_substrate_carbons``         - carbon count of the detected ligand.
* ``substrate_rgyr``              - radius of gyration (A) of the ligand heavy atoms (compactness).
* ``foldback_c1_to_distal``       - min distance (A) from C1 to a chain carbon >= FARCHAIN bonds
  away (small => chain curled toward C1 = cyclization-compatible); NaN if no distal carbon.
* ``substrate_endtoend``          - max distance (A) from C1 to any ligand carbon (extension).
* ``n_aromatic_carbon_contacts``  - ligand carbons within AROMATIC_CUTOFF of an aromatic centroid.
* ``frac_aromatic_track``         - that count / n_substrate_carbons (track coverage 0-1).
* ``n_aromatics_lining``          - distinct aromatic residues within AROMATIC_CUTOFF of any carbon.
* ``mean_carbon_to_aromatic``     - mean over ligand carbons of distance to the nearest aromatic.
* ``n_residues``                  - modelled residue count, for context.

Citation: Christianson, D. W. Chem. Rev. 2017, 117, 11570-11648 (class-I TPS carbocation
cyclization; aromatic carbocation stabilization).
"""

import argparse
import os
import warnings
from collections import deque
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from active_site_geometry import structure_sequence_residues_atoms
from substrate_positioning import read_substrate_ligand
from plddt import _collect_structures

# Aromatic ring atoms whose centroid defines the cation-pi face (Biopython names).
AROMATIC_RING_ATOMS: Dict[str, Tuple[str, ...]] = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
# Substrate-carbon -> aromatic-ring-centroid distance (A) counted as a cation-pi contact.
DEFAULT_AROMATIC_CUTOFF = 6.0
# A chain carbon this many C-C bonds (or more) from C1 counts as 'distal' for the fold-back test.
DEFAULT_FARCHAIN_BONDS = 6

COLUMNS = [
    "ID",
    "substrate_present",
    "n_substrate_carbons",
    "substrate_rgyr",
    "foldback_c1_to_distal",
    "substrate_endtoend",
    "n_aromatic_carbon_contacts",
    "frac_aromatic_track",
    "n_aromatics_lining",
    "mean_carbon_to_aromatic",
    "n_residues",
]


def _aromatic_centroids(residues: list) -> np.ndarray:
    """Ring-centroid (n, 3) for each aromatic protein residue (>=4 ring atoms present)."""
    cents: List[np.ndarray] = []
    for res in residues:
        names = AROMATIC_RING_ATOMS.get(res.get_resname().strip().upper())
        if not names:
            continue
        pts = [np.asarray(res[a].get_coord(), float) for a in names if a in res]
        if len(pts) >= 4:
            cents.append(np.mean(pts, axis=0))
    return np.vstack(cents) if cents else np.empty((0, 3))


def _chain_topological_distance(carbons: np.ndarray) -> np.ndarray:
    """BFS bond-distance from C1 (== index 0) to each carbon along C-C bonds (<1.8 A);
    unreachable carbons get -1."""
    n = len(carbons)
    d = np.sqrt(((carbons[:, None, :] - carbons[None, :, :]) ** 2).sum(2))
    adjacency = [np.where((d[i] < 1.8) & (d[i] > 0))[0] for i in range(n)]
    topo = [-1] * n
    topo[0] = 0
    queue = deque([0])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            if topo[v] < 0:
                topo[v] = topo[u] + 1
                queue.append(v)
    return np.array(topo)


def cyclization_geometry(
    structure_path: str,
    *,
    aromatic_cutoff: float = DEFAULT_AROMATIC_CUTOFF,
    farchain_bonds: int = DEFAULT_FARCHAIN_BONDS,
    ion_resnames: Tuple[str, ...] = ("MG", "MN"),
    min_substrate_carbons: int = 5,
    substrate_resname: Optional[str] = None,
) -> Dict[str, object]:
    """Cyclization-relevant geometry for one structure. No prenyl-PP substrate -> a graceful
    not-applicable row (substrate_present=False, geometry NaN)."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, object] = {
        "substrate_present": False,
        "n_substrate_carbons": 0,
        "substrate_rgyr": np.nan,
        "foldback_c1_to_distal": np.nan,
        "substrate_endtoend": np.nan,
        "n_aromatic_carbon_contacts": 0,
        "frac_aromatic_track": np.nan,
        "n_aromatics_lining": 0,
        "mean_carbon_to_aromatic": np.nan,
        "n_residues": len(sequence),
    }

    resname, lig_coords, _lig_b, diphos, carbons, _ions = read_substrate_ligand(
        structure_path, ion_resnames=ion_resnames,
        min_carbons=min_substrate_carbons, substrate_resname=substrate_resname,
    )
    if not len(carbons) or not len(diphos):
        return result
    result["substrate_present"] = True
    result["n_substrate_carbons"] = int(len(carbons))

    # C1 = ligand carbon nearest the diphosphate (same definition as substrate_positioning).
    d_c_to_diphos = np.sqrt(((carbons[:, None, :] - diphos[None, :, :]) ** 2).sum(2)).min(1)
    carbons = carbons[np.argsort(d_c_to_diphos)]   # C1 first
    c1 = carbons[0]

    result["substrate_rgyr"] = float(np.sqrt(((lig_coords - lig_coords.mean(0)) ** 2).sum(1).mean()))
    topo = _chain_topological_distance(carbons)
    far = np.where(topo >= farchain_bonds)[0]
    if len(far):
        result["foldback_c1_to_distal"] = float(np.sqrt(((carbons[far] - c1) ** 2).sum(1)).min())
    result["substrate_endtoend"] = float(np.sqrt(((carbons - c1) ** 2).sum(1)).max())

    aro = _aromatic_centroids(residues)
    if len(aro):
        dca = np.sqrt(((carbons[:, None, :] - aro[None, :, :]) ** 2).sum(2))   # carbons x aromatics
        nearest = dca.min(1)
        n_contact = int((nearest <= aromatic_cutoff).sum())
        result["n_aromatic_carbon_contacts"] = n_contact
        result["frac_aromatic_track"] = float(n_contact / len(carbons))
        result["mean_carbon_to_aromatic"] = float(nearest.mean())
        result["n_aromatics_lining"] = int((dca.min(0) <= aromatic_cutoff).sum())
    return result


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_cyclization_geometry.csv")


def cyclization_geometry_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    aromatic_cutoff: float = DEFAULT_AROMATIC_CUTOFF,
    farchain_bonds: int = DEFAULT_FARCHAIN_BONDS,
    ion_resnames: Tuple[str, ...] = ("MG", "MN"),
    min_substrate_carbons: int = 5,
    substrate_resname: Optional[str] = None,
) -> pd.DataFrame:
    """Cyclization-relevant geometry for every structure in structs_dir; CSV keyed by ID.
    Auto-detects an AF3 af_output dir vs a flat .pdb/.cif dir (plddt._collect_structures);
    writes <structs_dir>_cyclization_geometry.csv by default. No-substrate -> NaN row."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output dir with "
            "<job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif).")
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")
    print(f"Aromatic cutoff {aromatic_cutoff} A; far-chain >= {farchain_bonds} bonds; "
          f"min substrate carbons {min_substrate_carbons}"
          + (f"; forced substrate resname {substrate_resname}" if substrate_resname else ""))

    rows: List[Dict[str, object]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = cyclization_geometry(
                path, aromatic_cutoff=aromatic_cutoff, farchain_bonds=farchain_bonds,
                ion_resnames=ion_resnames, min_substrate_carbons=min_substrate_carbons,
                substrate_resname=substrate_resname)
        except Exception as exc:
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {c: (False if c == "substrate_present"
                         else 0 if c in ("n_substrate_carbons", "n_aromatic_carbon_contacts",
                                         "n_aromatics_lining", "n_residues")
                         else np.nan) for c in COLUMNS if c != "ID"}
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)
    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    n_sub = int(df["substrate_present"].sum())
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    print(f"  [info] {n_sub}/{len(df)} structures carry a prenyl-PP substrate (0 => apo / Mg-only / "
          "mg_ppi: cyclization_geometry is not applicable).")
    return df
