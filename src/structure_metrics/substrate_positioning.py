"""Substrate-positioning (catalytic-site) check for AF3 holo-folded class-I TPS designs.

When AF3 co-folds a design WITH its prenyl-diphosphate substrate (``--af3_cofold mg_<sub>``
forces one substrate for all designs; ``mg_ee`` uses each design's EnzymeExplorer-predicted
substrate), this tool asks whether that substrate is actually POISED FOR CATALYSIS in the
carboxylate cage:

  * its DIPHOSPHATE (the leaving group) should sit at the DDXXD/NSE metal cage (bridged by the
    Mg cluster), and
  * its REACTIVE CARBON (C1 — where the allylic carbocation forms once PPi departs) should be
    held near that machinery, not flung into bulk solvent.

Like ``ion_site_check`` (and unlike the apo-robust active-site tools), it READS the ligand
HETATMs. The substrate ligand is AUTO-DETECTED by composition — the HETATM residue with >=1
phosphorus AND >= ``min_substrate_carbons`` carbons (a prenyl-PP) — which distinguishes it from
the bare ``POP`` pyrophosphate of ``mg_ppi`` (no carbons) and the monatomic Mg/Mn ions. For a
linear prenyl-PP (GPP/FPP/GGPP/GFPP) every oxygen is in the diphosphate, so the diphosphate
atom set is simply the ligand's P + O atoms; the reactive carbon is the ligand carbon nearest
the diphosphate. ``--substrate_resname`` forces a specific ligand residue name instead.

Columns (keyed by ``ID``):
* ``metal_point_found``               — bool; the DDXXD(+NSE/DTE) cage centroid is computable.
* ``substrate_present``               — bool; a prenyl-PP substrate ligand was found.
* ``substrate_resname``               — residue name of the detected ligand ('' if none).
* ``n_substrate_atoms``               — heavy-atom count of the ligand.
* ``substrate_plddt``                 — mean B-factor (pLDDT) over the ligand atoms (placement
  confidence); NaN if no ligand.
* ``diphosphate_to_cage_dist``        — diphosphate centroid -> cage centroid (A); NaN if no
  substrate or no metal point.
* ``min_diphosphate_to_cage_oxygen``  — closest diphosphate atom -> any cage carboxylate O (A).
* ``diphosphate_to_nearest_ion``      — closest diphosphate atom -> nearest Mg/Mn ion (A); NaN
  if no ions modelled.
* ``reactive_carbon_to_cage_dist``    — reactive carbon (C1) -> cage centroid (A).
* ``substrate_in_site`` (bool)        — diphosphate centroid within ``site_radius`` of the cage.
* ``n_residues``                      — modelled residue count, for context.

SEMANTICS / NOT-APPLICABLE: apo / ESMFold / Mg-only / mg_ppi structures carry no prenyl-PP
substrate, so ``substrate_present=False`` and the geometry columns are NaN — a graceful
not-applicable row. The CSV is still written in full and the tool exits 0. This tool therefore
only carries signal for SUBSTRATE holo folds (``--af3_cofold mg_<sub>`` / ``mg_ee``).

Citation: Christianson, D. W. Chem. Rev. 2017, 117, 11570-11648 (class-I TPS catalysis:
metal-triggered ionization of the prenyl-diphosphate substrate).
"""
from __future__ import annotations

import argparse
import os
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from active_site_geometry import (
    _coordinating_oxygens,
    coordinating_indices_relaxed,
    metal_point as _cage_metal_point,
    structure_sequence_residues_atoms,
    _parser_for,
)
from plddt import _collect_structures

# Ions to exclude when auto-detecting the substrate, and to measure the diphosphate->ion
# distance against.
DEFAULT_ION_RESNAMES = ("MG", "MN")
# Metal-ion ELEMENTS. Ion detection is element-based (not just resname) so it is robust to
# naming across structure sources: AF3 names the ions by resname 'MG', Boltz2 names them
# 'LIG2' (but the element column is still MG). A monatomic residue whose atoms are all these
# elements is treated as an ion regardless of its resname.
ION_ELEMENTS = {"MG", "MN"}
# A HETATM residue is a candidate prenyl-PP substrate if it has >=1 P and >= this many C.
DEFAULT_MIN_SUBSTRATE_CARBONS = 5
# Diphosphate centroid within this distance (A) of the cage centroid -> reported as
# diphosphate_to_cage_dist (informational; inherits the shared metal_point's oxygen-centroid,
# which can be pulled off-site when the relaxed coordinating set spans a splayed cage).
DEFAULT_SITE_RADIUS = 6.0
# The ROBUST "is the diphosphate at the cage" test: closest diphosphate atom -> any cage
# carboxylate oxygen within this distance (A). Direct PPi...Asp(-Mg) contact is ~2.5-3.5 A;
# 4.0 is a slightly relaxed cutoff. substrate_in_site is based on THIS, not the centroid.
DEFAULT_COORD_CUTOFF = 4.0
# Non-substrate HETATMs to always ignore.
_IGNORE_RESNAMES = {"HOH", "WAT", "POP", "PPV", "PPK"}

COLUMNS = [
    "ID",
    "metal_point_found",
    "substrate_present",
    "substrate_resname",
    "n_substrate_atoms",
    "substrate_plddt",
    "diphosphate_to_cage_dist",
    "min_diphosphate_to_cage_oxygen",
    "diphosphate_to_nearest_ion",
    "diphosphate_to_ion_centroid",
    "reactive_carbon_to_nearest_ion",
    "reactive_carbon_to_ion_centroid",
    "reactive_carbon_to_cage_dist",
    "substrate_in_site",
    "n_residues",
]


def _element(atom) -> str:
    el = (atom.element or "").strip().upper()
    if el:
        return el
    # Fallback: first alphabetic char(s) of the atom name (e.g. "C1" -> "C").
    name = atom.get_name().strip()
    return "".join(c for c in name if c.isalpha())[:1].upper()


def read_substrate_ligand(
    structure_path: str,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    min_carbons: int = DEFAULT_MIN_SUBSTRATE_CARBONS,
    substrate_resname: Optional[str] = None,
) -> Tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find the prenyl-PP substrate ligand among the HETATMs.

    Returns ``(resname, all_coords, all_bfactors, diphosphate_coords, carbon_coords,
    ion_coords)``. The substrate is the HETATM residue with the MOST carbons among residues
    having >=1 P and >= ``min_carbons`` C (or the residue named ``substrate_resname`` if given).
    Empty arrays / '' resname when no substrate is found. ``diphosphate_coords`` = the
    substrate's P + O atoms; ``carbon_coords`` = its C atoms; ``ion_coords`` = all Mg/Mn ion
    atoms (for the diphosphate->ion distance)."""
    ion_set = {r.strip().upper() for r in ion_resnames}
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))

    ions: List[np.ndarray] = []
    candidates: List[Tuple[str, List[Tuple[str, np.ndarray, float]]]] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":               # protein residue -> skip
                continue
            resname = residue.get_resname().strip().upper()
            atoms = [(_element(a), np.asarray(a.get_coord(), float), float(a.get_bfactor()))
                     for a in residue]
            if (resname in ion_set) or (atoms and all(el in ION_ELEMENTS for (el, _c, _b) in atoms)):
                ions.extend(c for (_el, c, _b) in atoms)  # element-based: AF3 'MG' & Boltz2 'LIG2'
                continue
            if resname in _IGNORE_RESNAMES:
                continue
            if substrate_resname is not None:
                if resname == substrate_resname.strip().upper():
                    candidates.append((resname, atoms))
                continue
            n_p = sum(el == "P" for (el, _c, _b) in atoms)
            n_c = sum(el == "C" for (el, _c, _b) in atoms)
            if n_p >= 1 and n_c >= min_carbons:
                candidates.append((resname, atoms))

    ion_arr = np.vstack(ions) if ions else np.empty((0, 3))
    empty = np.empty((0, 3))
    if not candidates:
        return "", empty, np.empty((0,)), empty, empty, ion_arr
    # Most carbons = the prenyl-PP substrate (over e.g. a small cofactor).
    resname, atoms = max(candidates, key=lambda rc: sum(el == "C" for (el, _c, _b) in rc[1]))
    all_coords = np.vstack([c for (_el, c, _b) in atoms])
    all_b = np.array([b for (_el, _c, b) in atoms], float)
    diphos = np.vstack([c for (el, c, _b) in atoms if el in ("P", "O")]) \
        if any(el in ("P", "O") for (el, _c, _b) in atoms) else empty
    carbons = np.vstack([c for (el, c, _b) in atoms if el == "C"]) \
        if any(el == "C" for (el, _c, _b) in atoms) else empty
    return resname, all_coords, all_b, diphos, carbons, ion_arr


def substrate_positioning(
    structure_path: str,
    *,
    site_radius: float = DEFAULT_SITE_RADIUS,
    coord_cutoff: float = DEFAULT_COORD_CUTOFF,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    min_substrate_carbons: int = DEFAULT_MIN_SUBSTRATE_CARBONS,
    substrate_resname: Optional[str] = None,
) -> Dict[str, object]:
    """Substrate-positioning metrics for one structure. No substrate ligand -> a graceful
    not-applicable row (substrate_present=False, geometry NaN)."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, object] = {
        "metal_point_found": False,
        "substrate_present": False,
        "substrate_resname": "",
        "n_substrate_atoms": 0,
        "substrate_plddt": np.nan,
        "diphosphate_to_cage_dist": np.nan,
        "min_diphosphate_to_cage_oxygen": np.nan,
        "diphosphate_to_nearest_ion": np.nan,
        "diphosphate_to_ion_centroid": np.nan,
        "reactive_carbon_to_nearest_ion": np.nan,
        "reactive_carbon_to_ion_centroid": np.nan,
        "reactive_carbon_to_cage_dist": np.nan,
        "substrate_in_site": False,
        "n_residues": len(sequence),
    }

    resname, lig_coords, lig_b, diphos, carbons, ions = read_substrate_ligand(
        structure_path, ion_resnames=ion_resnames,
        min_carbons=min_substrate_carbons, substrate_resname=substrate_resname,
    )
    if len(lig_coords):
        result["substrate_present"] = True
        result["substrate_resname"] = resname
        result["n_substrate_atoms"] = int(len(lig_coords))
        result["substrate_plddt"] = float(np.mean(lig_b)) if len(lig_b) else np.nan
        if len(diphos) and len(ions):
            result["diphosphate_to_nearest_ion"] = float(
                np.sqrt(((diphos[:, None, :] - ions[None, :, :]) ** 2).sum(2)).min())

    # Ion-anchored substrate geometry (reference-independent: measured against the cofolded
    # Mg cluster, NOT the apo metal_point which mislocalizes for two-domain folds). Reactive
    # carbon = ligand carbon nearest the diphosphate (same definition as the cage version).
    if len(ions):
        ion_centroid = ions.mean(axis=0)
        if len(diphos):
            result["diphosphate_to_ion_centroid"] = float(np.sqrt(((diphos.mean(axis=0) - ion_centroid) ** 2).sum()))
        if len(carbons) and len(diphos):
            _d_c = np.sqrt(((carbons[:, None, :] - diphos[None, :, :]) ** 2).sum(2)).min(1)
            _reactive_c = carbons[int(np.argmin(_d_c))]
            result["reactive_carbon_to_nearest_ion"] = float(np.sqrt(((_reactive_c[None, :] - ions) ** 2).sum(1)).min())
            result["reactive_carbon_to_ion_centroid"] = float(np.sqrt(((_reactive_c - ion_centroid) ** 2).sum()))

    cage = _cage_metal_point(sequence, residues)
    if cage is None or not len(diphos):
        return result                              # no anchor or no substrate -> NaN geometry
    result["metal_point_found"] = True

    diphos_centroid = diphos.mean(axis=0)
    result["diphosphate_to_cage_dist"] = float(np.sqrt(((diphos_centroid - cage) ** 2).sum()))

    idx = coordinating_indices_relaxed(sequence)
    cage_oxygens = _coordinating_oxygens(idx, residues) if idx is not None else np.empty((0, 3))
    if len(cage_oxygens):
        min_o = float(np.sqrt(((diphos[:, None, :] - cage_oxygens[None, :, :]) ** 2).sum(2)).min())
        result["min_diphosphate_to_cage_oxygen"] = min_o
        # Robust in-site test: the diphosphate reaches a cage carboxylate oxygen. (The
        # centroid distance above can be inflated by a splayed coordinating-oxygen set.)
        result["substrate_in_site"] = bool(min_o <= coord_cutoff)

    # Reactive carbon = ligand carbon nearest the diphosphate (the C bonded to the bridging
    # ester O); measure its distance to the cage centroid.
    if len(carbons):
        d_c_to_diphos = np.sqrt(((carbons[:, None, :] - diphos[None, :, :]) ** 2).sum(2)).min(1)
        reactive_c = carbons[int(np.argmin(d_c_to_diphos))]
        result["reactive_carbon_to_cage_dist"] = float(np.sqrt(((reactive_c - cage) ** 2).sum()))
    return result


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_substrate_positioning.csv")


def substrate_positioning_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    site_radius: float = DEFAULT_SITE_RADIUS,
    coord_cutoff: float = DEFAULT_COORD_CUTOFF,
    ion_resnames: Tuple[str, ...] = DEFAULT_ION_RESNAMES,
    min_substrate_carbons: int = DEFAULT_MIN_SUBSTRATE_CARBONS,
    substrate_resname: Optional[str] = None,
) -> pd.DataFrame:
    """Substrate-positioning check for every structure in ``structs_dir``; CSV keyed by ID.
    Auto-detects an AF3 ``af_output`` dir vs a flat dir of .pdb/.cif (via
    ``plddt._collect_structures``); writes ``<structs_dir>_substrate_positioning.csv`` by
    default. Structures with no prenyl-PP substrate get a graceful not-applicable row."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output dir with "
            "<job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif).")
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")
    print(f"Site radius {site_radius} A; ions {sorted({r.upper() for r in ion_resnames})}; "
          f"min substrate carbons {min_substrate_carbons}"
          + (f"; forced substrate resname {substrate_resname}" if substrate_resname else ""))

    rows: List[Dict[str, object]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = substrate_positioning(
                path, site_radius=site_radius, coord_cutoff=coord_cutoff, ion_resnames=ion_resnames,
                min_substrate_carbons=min_substrate_carbons, substrate_resname=substrate_resname)
        except Exception as exc:
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {c: (False if c in ("metal_point_found", "substrate_present", "substrate_in_site")
                         else "" if c == "substrate_resname"
                         else 0 if c in ("n_substrate_atoms", "n_residues") else np.nan)
                     for c in COLUMNS if c != "ID"}
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
    print(f"  [info] {n_sub}/{len(df)} structures carry a prenyl-PP substrate ligand "
          "(0 => apo / Mg-only / mg_ppi: substrate_positioning is not applicable).")
    return df
