from __future__ import annotations

"""Aromatic / cation-pi pocket-lining metric for class-I TPS designs.

Carbocation intermediates along a terpene-cyclization cascade are stabilized by
the quadrupole of aromatic side chains (Trp >> Tyr ~ Phe) that line the catalytic
pocket — the classic "aromatic box". The *count and orientation* of those pocket
aromatics is therefore a fold-agnostic proxy for cyclization capability: a flat,
aromatic-poor pocket can still bind the prenyl-diphosphate substrate but is far
less able to template a polycyclic product.

This tool reports RAW numbers only (no "natural TPS" band — that's the reference-
stats pipeline's job). It is fold-agnostic (AlphaFold / ESMFold .pdb/.cif) and
apo-robust (no metals / substrate need be modelled), but see the apo caveat below.

For each structure in ``structs_dir`` we

1. parse it with Biopython (REUSING the af_output/flat auto-detection, the ID-stem
   convention, and the sequence/residue/atom derivation from the neighbouring
   structure tools — ``plddt._collect_structures`` / ``active_site_geometry.
   structure_sequence_residues_atoms``),
2. derive the **carboxylate-cage metal point** via the canonical, single-source-of-
   truth ``active_site_geometry.metal_point`` = centroid of the DDXXD (+ NSE/DTE when
   that motif is also matched) coordinating side-chain oxygens. DDXXD is required;
   NSE/DTE is an optional refinement, so real DDXXD-only TPS (e.g. TEAS/5EAT) get a
   metal point instead of an all-NaN row,
3. select the **pocket residues** = every residue with ANY atom within ``--cutoff``
   (default 10 A) of that metal point (a distance shell — fold-agnostic, no
   pocket-detection geometry assumed), and
4. among those pocket residues compute the metrics below.

POCKET-SELECTION APPROACH
-------------------------
A spherical distance shell around the carboxylate-cage metal point. The metal point
is where the trinuclear Mg/Mn cluster + the diphosphate of the substrate sit, so the
prenyl chain (and the migrating carbocation) extends into the hydrophobic cavity just
"below" it. A ~10 A shell captures the residues whose side chains line that cavity
mouth and walls. The cutoff is a CLI arg so the reference-stats pipeline can sweep it.

RING-ORIENTATION APPROACH (``n_inward_facing_aromatics``)
---------------------------------------------------------
For each pocket aromatic we build the six-membered ring (Phe/Tyr: CG CD1 CD2 CE1 CE2
CZ; Trp uses its six-membered benzene ring CD2 CE2 CE3 CZ2 CZ3 CH2), take the ring
centroid and the ring-plane normal (SVD of the centred ring atoms — smallest singular
vector). A ring "faces" the cavity interior when (a) its centroid is within cation-pi
range of the cation locus AND (b) the ring FACE (not the edge) points toward the locus,
i.e. the angle between the ring normal and the centroid->locus vector is small (the pi
face, perpendicular to the normal, is what stacks with the cation). Defaults: distance
3.5-6.0 A (``--cation_pi_min`` / ``--cation_pi_max``); face-on angle <= 45 deg from the
normal (``--face_angle_deg``); normal sign is folded (a ring has two equivalent faces),
so we use ``abs(cos)``.

CATION-LOCUS APPROXIMATION
--------------------------
The true carbocation locus is unknown without the bound intermediate. We approximate
it as the **pocket-residue Calpha centroid** (the geometric centre of the residues
lining the pocket) — a stable, apo-available stand-in for "the middle of the cavity".
This is offset from, and complementary to, the metal point (which sits at the cavity
mouth where the diphosphate anchors): the carbocation migrates DOWN into the cavity
body, which the Calpha centroid of the lining residues tracks. Documented as an
approximation; ``n_pocket_aromatics`` (the count) does not depend on it.

APO CAVEAT
----------
Counts (``n_pocket_aromatics`` / per-type / ``aromatic_fraction``) are robust: they
depend only on residue identity + Calpha-ish proximity, which a confident apo predicted
model gets right. ``n_inward_facing_aromatics`` is ROTAMER-SENSITIVE: side-chain ring
orientation in an apo model is a prediction, not a fact, so treat the inward-facing
count as a softer signal than the raw aromatic count.

Metrics (columns, keyed by ``ID``):

* ``metal_point_found`` — bool; False (and counts NaN) when DDXXD or its coordinating
  oxygens are missing so no metal point could be placed. A real RED FLAG for a class-I
  design — recorded, never silently dropped. Every structure gets exactly one row.
* ``n_pocket_residues`` — residues with any atom within ``cutoff`` of the metal point.
* ``n_pocket_aromatics`` — Trp + Tyr + Phe among the pocket residues (core count).
* ``n_trp`` / ``n_tyr`` / ``n_phe`` — per-type breakdown (Trp is the strongest
  cation-pi stabilizer; reported separately).
* ``n_his`` — His among the pocket residues (also cation-pi-capable; kept OUT of the
  core aromatic count, reported separately).
* ``aromatic_fraction`` — ``n_pocket_aromatics / n_pocket_residues``.
* ``n_inward_facing_aromatics`` — pocket aromatics whose ring face points at the
  cation locus and lies within cation-pi range (the geometry refinement above).
* ``n_residues`` — modelled residue count, for context.
"""

import argparse
import os
import sys
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Reuse the shared structure loader + the carboxylate-cage machinery (single source
# of truth for the metal point) from the neighbouring structure tools.
from active_site_geometry import (
    metal_point as _cage_metal_point,
    structure_sequence_residues_atoms,
)
from plddt import _collect_structures

# Defaults (CLI-overridable). Cutoff is the pocket-shell radius; the cation-pi
# window + face angle govern the inward-facing geometry refinement.
DEFAULT_CUTOFF = 10.0
DEFAULT_CATION_PI_MIN = 3.5
DEFAULT_CATION_PI_MAX = 6.0
DEFAULT_FACE_ANGLE_DEG = 45.0

# Core cation-pi aromatic residues (Trp is reported separately as the strongest).
AROMATIC_RESNAMES = ("TRP", "TYR", "PHE")

# Six-membered aromatic ring atoms used for centroid + plane normal, per residue.
# Phe/Tyr: the benzene ring. Trp: the six-membered benzene ring of the indole
# (the larger, more polarizable face that dominates cation-pi stacking).
RING_ATOMS: Dict[str, Tuple[str, ...]] = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
}

COLUMNS = [
    "ID",
    "metal_point_found",
    "n_pocket_residues",
    "n_pocket_aromatics",
    "n_trp",
    "n_tyr",
    "n_phe",
    "n_his",
    "aromatic_fraction",
    "n_inward_facing_aromatics",
    "n_residues",
]


def _residue_coords(residue) -> np.ndarray:
    """(n_atoms, 3) coordinates of one residue's atoms."""
    return np.vstack([np.asarray(a.get_coord(), dtype=float) for a in residue])


def _pocket_residue_indices(
    metal_point: np.ndarray, residues: list, cutoff: float
) -> List[int]:
    """Indices of residues with ANY atom within ``cutoff`` of the metal point."""
    pocket: List[int] = []
    cutoff_sq = cutoff * cutoff
    for i, residue in enumerate(residues):
        coords = _residue_coords(residue)
        if coords.size == 0:
            continue
        if (((coords - metal_point) ** 2).sum(axis=1) <= cutoff_sq).any():
            pocket.append(i)
    return pocket


def _ring_centroid_normal(residue) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """(centroid, unit normal) of an aromatic residue's six-membered ring, or None
    if the ring atoms aren't all present. Normal = smallest right-singular vector of
    the centred ring coordinates (the ring-plane normal)."""
    names = RING_ATOMS.get(residue.get_resname())
    if names is None:
        return None
    pts: List[np.ndarray] = []
    for name in names:
        if name not in residue:
            return None
        pts.append(np.asarray(residue[name].get_coord(), dtype=float))
    ring = np.vstack(pts)
    centroid = ring.mean(axis=0)
    _, _, vh = np.linalg.svd(ring - centroid)
    normal = vh[-1]
    norm = np.linalg.norm(normal)
    if norm == 0:
        return None
    return centroid, normal / norm


def _cation_locus(pocket_indices: List[int], residues: list) -> Optional[np.ndarray]:
    """Carbocation-locus approximation = Calpha centroid of the pocket residues
    (the geometric centre of the lining residues). Falls back to CB / any atom when
    a residue lacks CA. None if no usable point exists."""
    pts: List[np.ndarray] = []
    for i in pocket_indices:
        residue = residues[i]
        if "CA" in residue:
            pts.append(np.asarray(residue["CA"].get_coord(), dtype=float))
        elif "CB" in residue:
            pts.append(np.asarray(residue["CB"].get_coord(), dtype=float))
    if not pts:
        return None
    return np.vstack(pts).mean(axis=0)


def aromatic_lining(
    structure_path: str,
    *,
    cutoff: float = DEFAULT_CUTOFF,
    cation_pi_min: float = DEFAULT_CATION_PI_MIN,
    cation_pi_max: float = DEFAULT_CATION_PI_MAX,
    face_angle_deg: float = DEFAULT_FACE_ANGLE_DEG,
) -> Dict[str, float]:
    """Aromatic / cation-pi pocket-lining metrics for one structure. Counts are NaN
    and ``metal_point_found`` False when the metal point can't be placed."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, float] = {
        "metal_point_found": False,
        "n_pocket_residues": np.nan,
        "n_pocket_aromatics": np.nan,
        "n_trp": np.nan,
        "n_tyr": np.nan,
        "n_phe": np.nan,
        "n_his": np.nan,
        "aromatic_fraction": np.nan,
        "n_inward_facing_aromatics": np.nan,
        "n_residues": len(sequence),
    }

    metal_point = _cage_metal_point(sequence, residues)
    if metal_point is None:
        return result
    result["metal_point_found"] = True

    pocket_idx = _pocket_residue_indices(metal_point, residues, cutoff)
    n_pocket = len(pocket_idx)
    result["n_pocket_residues"] = n_pocket

    n_trp = sum(residues[i].get_resname() == "TRP" for i in pocket_idx)
    n_tyr = sum(residues[i].get_resname() == "TYR" for i in pocket_idx)
    n_phe = sum(residues[i].get_resname() == "PHE" for i in pocket_idx)
    n_his = sum(residues[i].get_resname() == "HIS" for i in pocket_idx)
    n_aromatics = n_trp + n_tyr + n_phe
    result["n_trp"] = n_trp
    result["n_tyr"] = n_tyr
    result["n_phe"] = n_phe
    result["n_his"] = n_his
    result["n_pocket_aromatics"] = n_aromatics
    result["aromatic_fraction"] = (n_aromatics / n_pocket) if n_pocket else np.nan

    # Inward-facing refinement against the cation-locus approximation.
    locus = _cation_locus(pocket_idx, residues)
    n_inward = 0
    if locus is not None:
        cos_thresh = np.cos(np.radians(face_angle_deg))
        for i in pocket_idx:
            residue = residues[i]
            if residue.get_resname() not in AROMATIC_RESNAMES:
                continue
            cn = _ring_centroid_normal(residue)
            if cn is None:
                continue
            centroid, normal = cn
            to_locus = locus - centroid
            dist = float(np.linalg.norm(to_locus))
            if not (cation_pi_min <= dist <= cation_pi_max):
                continue
            # Face-on: the centroid->locus vector aligns with the ring normal
            # (the pi face is perpendicular to the ring plane). Fold the normal
            # sign (a ring presents two equivalent faces).
            cos_face = abs(float(np.dot(normal, to_locus) / dist))
            if cos_face >= cos_thresh:
                n_inward += 1
    result["n_inward_facing_aromatics"] = n_inward
    return result


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_aromatic_lining.csv")


def aromatic_lining_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    cutoff: float = DEFAULT_CUTOFF,
    cation_pi_min: float = DEFAULT_CATION_PI_MIN,
    cation_pi_max: float = DEFAULT_CATION_PI_MAX,
    face_angle_deg: float = DEFAULT_FACE_ANGLE_DEG,
) -> pd.DataFrame:
    """Aromatic / cation-pi pocket-lining for every structure in ``structs_dir``;
    CSV keyed by ID (filename stem / af3 job name). Every structure gets one row;
    structures with no locatable metal point are kept with ``metal_point_found``
    False and NaN counts (a recorded red flag)."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")
    print(
        f"Pocket cutoff {cutoff} A; cation-pi window {cation_pi_min}-{cation_pi_max} A; "
        f"face angle <= {face_angle_deg} deg."
    )

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = aromatic_lining(
                path,
                cutoff=cutoff,
                cation_pi_min=cation_pi_min,
                cation_pi_max=cation_pi_max,
                face_angle_deg=face_angle_deg,
            )
        except Exception as exc:  # malformed/unparsable -> red-flag row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {
                "metal_point_found": False,
                "n_pocket_residues": np.nan,
                "n_pocket_aromatics": np.nan,
                "n_trp": np.nan,
                "n_tyr": np.nan,
                "n_phe": np.nan,
                "n_his": np.nan,
                "aromatic_fraction": np.nan,
                "n_inward_facing_aromatics": np.nan,
                "n_residues": 0,
            }
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    n_no_metal = int((~df["metal_point_found"]).sum())
    if n_no_metal:
        print(
            f"  [note] {n_no_metal}/{len(df)} structure(s) had no locatable metal point "
            "(missing DDXXD motif or its coordinating oxygens) -> counts NaN."
        )

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    return df
