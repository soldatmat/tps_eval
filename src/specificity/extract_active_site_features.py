from __future__ import annotations

"""Active-site / cation-specific-residue feature vectors for MARTS-DB TPS.

Goal
----
Build a FIXED-LENGTH, per-protein feature vector describing the catalytic
active-site region of each terpene synthase, intended for a downstream UMAP map
of first-cyclization product-class structure. The premise (Durairaj et al.,
PLOS Comp Biol 2021): TPS first-cyclization specificity is dominated by the
active-site contour residues lining the Mg2+-binding carboxylate cage (the
DDXXD + NSE/DTE motif neighbourhood), NOT the global fold. A UMAP on
"cation-specific" active-site residues separates sesquiterpene product classes
that a whole-sequence/whole-fold embedding does not.

Why a property-PROFILE vector (and not an aligned-position one-hot)
-------------------------------------------------------------------
The reference set spans 22 product classes across mono/sesqui/di/sester-terpene
class-I synthases at <20% global identity and multiple folds. Building a
*consistent, aligned* set of active-site positions across all of them (the
Durairaj approach uses a structural superposition + curated cation-specific
positions on a single fold family) is the crux and is NOT tractable to do
robustly here. Instead we sidestep alignment entirely:

1. Anchor on the **carboxylate-cage metal point** — the centroid of the DDXXD
   (+ NSE/DTE when matched) coordinating side-chain oxygens — via the CANONICAL
   single-source-of-truth ``active_site_geometry.metal_point`` (DDXXD required;
   NSE/DTE optional refinement; same definition used by ``aromatic_lining`` /
   ``pocket_descriptors``). This is the one physically-defined, fold-agnostic
   reference point shared by every class-I TPS active site.
2. Select the **active-site shell** = residues with any atom within ``--radius``
   (default 12 A) of the metal point. This captures the pocket walls / cavity
   mouth without assuming any alignment or pocket-detection geometry.
3. Featurize the shell as an ALIGNMENT-FREE, fixed-length descriptor: a
   physico-chemical PROPERTY PROFILE (fraction of shell residues in each of
   several property classes) plus a handful of geometric / cage descriptors.
   This is order-invariant and length-invariant, so it is comparable across
   folds and classes without a multiple alignment.

The OSC / triterpene outlier (product class 12)
-----------------------------------------------
Class-12 enzymes are oxidosqualene cyclases: a DIFFERENT fold, NO Mg2+ cluster,
NO DDXXD/NSE-DTE motifs (substrate is (S)-2,3-epoxysqualene, protonation-
initiated, not metal/diphosphate-ionization-initiated). They therefore get NO
locatable metal point and an all-NaN feature row, recorded with
``metal_point_found=False``. They are thus INTRINSICALLY excluded from the
active-site feature space — which is correct: there is no common active-site
frame between the class-I TPS fold and the OSC fold. The CSV keeps their rows
(flagged) so the downstream map can drop them explicitly.

Feature vector (columns)
-------------------------
Metadata / context (not features per se):
* ``id``                    — Enzyme_marts_ID (CSV key).
* ``metal_point_found``     — bool; False => OSC-outlier / no class-I cage =>
                              feature columns NaN.
* ``n_shell_residues``      — residues in the active-site shell.
* ``n_residues``            — modelled residue count.
* ``radius_A``              — the shell radius used (constant per run).

Composition PROPERTY PROFILE (fractions over the shell residues; sum of the
mutually-exclusive property partition is 1):
* ``frac_aromatic``         — F,Y,W (cation-pi stabilizers — the dominant TPS
                              active-site signal).
* ``frac_aliphatic``        — A,V,L,I,M (hydrophobic cavity walls).
* ``frac_acidic``           — D,E (the metal-binding carboxylates + acid/base).
* ``frac_basic``            — K,R,H (cation-binding / diphosphate-anchoring).
* ``frac_polar``            — S,T,N,Q,C (H-bonding / nucleophile-positioning).
* ``frac_glycine``          — G (backbone flexibility / cavity shaping).
* ``frac_proline``          — P (loop rigidity).

Per-amino-acid fractions (20 columns ``frac_aa_<X>``) — the full residue
composition of the shell, so the downstream map can learn finer distinctions
than the 7 property buckets. Together these are the "cation-specific residue
profile".

Geometry / cage descriptors (physical pocket shape, alignment-free):
* ``carboxylate_convergence_radius`` — RMS spread of the coordinating oxygens
  (cage tightness; reused from active_site_geometry semantics).
* ``n_coordinating_oxygens``         — count of cage oxygens found.
* ``shell_radius_of_gyration``       — Rg of the shell CA atoms about their
  centroid (active-site cavity scale; tracks substrate size GPP<FPP<GGPP<GFPP).
* ``mean_dist_to_metal_point``       — mean CA distance of shell residues to the
  metal point (pocket compactness).
* ``n_aromatic_within_8A``           — aromatics whose CA is within 8 A of the
  metal point (the tight cation-pi inner shell — a sharper specificity signal
  than the full 12 A count).

Dimensionality: 7 property fractions + 20 per-AA fractions + 5 geometric = 32
numeric feature dimensions (plus 5 metadata columns).

The feature columns are NaN for any protein with no locatable metal point
(class-12 OSC, or a degenerate/failed structure). RAW numbers only — no
class-conditional banding (that is the map's job).

CLI / dir conventions (af3-vs-flat detection, ID = filename stem) mirror
``plddt.py`` / ``active_site_geometry.py``. This module ONLY ADDS code; it
imports the existing motif/cage machinery read-only.
"""

import argparse
import os
import sys
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# Read-only reuse of the existing, agent-owned structure machinery.
from structure_metrics.active_site_geometry import (  # noqa: E402
    coordinating_indices_relaxed,
    metal_point as _cage_metal_point,
    structure_sequence_residues_atoms,
    _collect_structures,
    _coordinating_oxygens,
)

DEFAULT_RADIUS = 12.0
INNER_AROMATIC_RADIUS = 8.0

# Standard 20 amino acids (1-letter), fixed order for stable column ordering.
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# Mutually-exclusive physico-chemical partition of the 20 AAs.
PROPERTY_GROUPS: "OrderedDict[str, str]" = OrderedDict(
    [
        ("frac_aromatic", "FYW"),
        ("frac_aliphatic", "AVLIM"),
        ("frac_acidic", "DE"),
        ("frac_basic", "KRH"),
        ("frac_polar", "STNQC"),
        ("frac_glycine", "G"),
        ("frac_proline", "P"),
    ]
)

AROMATIC = set("FYW")

# Column layout.
META_COLUMNS = ["id", "metal_point_found", "n_shell_residues", "n_residues", "radius_A"]
PROFILE_COLUMNS = list(PROPERTY_GROUPS.keys())
AA_COLUMNS = [f"frac_aa_{a}" for a in AA_ORDER]
GEOM_COLUMNS = [
    "carboxylate_convergence_radius",
    "n_coordinating_oxygens",
    "shell_radius_of_gyration",
    "mean_dist_to_metal_point",
    "n_aromatic_within_8A",
]
FEATURE_COLUMNS = PROFILE_COLUMNS + AA_COLUMNS + GEOM_COLUMNS
COLUMNS = META_COLUMNS + FEATURE_COLUMNS

# Map Biopython 3-letter resname -> 1-letter. Reuse the same helper semantics as
# active_site_geometry (X for non-standard); but we already have the 1-letter
# sequence from structure_sequence_residues_atoms, so we index into that instead.


def _residue_ca(residue) -> Optional[np.ndarray]:
    if "CA" in residue:
        return np.asarray(residue["CA"].get_coord(), dtype=float)
    if "CB" in residue:
        return np.asarray(residue["CB"].get_coord(), dtype=float)
    return None


def _shell_indices(metal_point: np.ndarray, residues: list, radius: float) -> List[int]:
    """Indices of residues with ANY atom within ``radius`` of the metal point."""
    out: List[int] = []
    r2 = radius * radius
    for i, residue in enumerate(residues):
        coords = np.vstack(
            [np.asarray(a.get_coord(), dtype=float) for a in residue]
        ) if len(list(residue)) else np.empty((0, 3))
        if coords.size == 0:
            continue
        if (((coords - metal_point) ** 2).sum(axis=1) <= r2).any():
            out.append(i)
    return out


def _nan_features() -> Dict[str, float]:
    d: Dict[str, float] = {c: np.nan for c in FEATURE_COLUMNS}
    d["n_coordinating_oxygens"] = 0
    return d


def active_site_features(
    structure_path: str, *, radius: float = DEFAULT_RADIUS
) -> Dict[str, float]:
    """Active-site property/geometry feature vector for one structure.

    All feature columns are NaN (and ``metal_point_found`` False) when the
    carboxylate-cage metal point can't be placed (no DDXXD / no coordinating
    oxygen) — e.g. the OSC class-12 outlier."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, float] = {
        "metal_point_found": False,
        "n_shell_residues": 0,
        "n_residues": len(sequence),
        "radius_A": radius,
    }
    result.update(_nan_features())

    metal_point = _cage_metal_point(sequence, residues)
    if metal_point is None:
        return result
    result["metal_point_found"] = True

    shell = _shell_indices(metal_point, residues, radius)
    n_shell = len(shell)
    result["n_shell_residues"] = n_shell
    if n_shell == 0:
        return result

    # 1-letter codes of the shell residues (from the parsed sequence, aligned by
    # index to ``residues``). Non-standard residues are 'X' and are counted in
    # the denominator but contribute to no group/AA fraction.
    shell_aas = [sequence[i] for i in shell]

    # Per-AA composition fractions.
    counts = {a: 0 for a in AA_ORDER}
    for a in shell_aas:
        if a in counts:
            counts[a] += 1
    for a in AA_ORDER:
        result[f"frac_aa_{a}"] = counts[a] / n_shell

    # Property-group fractions (mutually exclusive partition of the standard 20).
    for col, members in PROPERTY_GROUPS.items():
        result[col] = sum(counts[a] for a in members) / n_shell

    # --- Geometry / cage descriptors ---
    # Carboxylate convergence radius + coordinating-oxygen count (cage tightness).
    cidx = coordinating_indices_relaxed(sequence)
    if cidx is not None:
        oxygens = _coordinating_oxygens(cidx, residues)
        result["n_coordinating_oxygens"] = int(len(oxygens))
        if len(oxygens):
            c = oxygens.mean(axis=0)
            result["carboxylate_convergence_radius"] = float(
                np.sqrt(((oxygens - c) ** 2).sum(axis=1).mean())
            )

    # Shell CA-based geometry.
    ca = [_residue_ca(residues[i]) for i in shell]
    ca = np.vstack([c for c in ca if c is not None]) if any(c is not None for c in ca) else np.empty((0, 3))
    if len(ca):
        centroid = ca.mean(axis=0)
        result["shell_radius_of_gyration"] = float(
            np.sqrt(((ca - centroid) ** 2).sum(axis=1).mean())
        )
        result["mean_dist_to_metal_point"] = float(
            np.sqrt(((ca - metal_point) ** 2).sum(axis=1)).mean()
        )

    # Inner-shell aromatics (CA within INNER_AROMATIC_RADIUS of the metal point).
    n_inner_arom = 0
    inner2 = INNER_AROMATIC_RADIUS * INNER_AROMATIC_RADIUS
    for i in shell:
        if sequence[i] in AROMATIC:
            ca_i = _residue_ca(residues[i])
            if ca_i is not None and ((ca_i - metal_point) ** 2).sum() <= inner2:
                n_inner_arom += 1
    result["n_aromatic_within_8A"] = n_inner_arom

    return result


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_active_site_features.csv")


def active_site_features_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    radius: float = DEFAULT_RADIUS,
    id_filter: Optional[set] = None,
) -> pd.DataFrame:
    """Active-site feature vectors for every structure in ``structs_dir``;
    CSV keyed by ``id`` (filename stem). ``id_filter`` (a set of IDs) restricts
    processing to those stems when given (e.g. the MARTS-DB enzyme set)."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    if id_filter is not None:
        structures = OrderedDict(
            (k, v) for k, v in structures.items() if k in id_filter
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) to process in {structs_dir}")
    print(f"Active-site shell radius: {radius} A")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = active_site_features(path, radius=radius)
        except Exception as exc:  # malformed/unparsable -> NaN row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {
                "metal_point_found": False,
                "n_shell_residues": 0,
                "n_residues": 0,
                "radius_A": radius,
            }
            stats.update(_nan_features())
            n_failed += 1
        stats["id"] = str(stem).strip()
        rows.append(stats)
        if i % 100 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("id").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    df.to_csv(save_path, index=False)

    n_no_metal = int((~df["metal_point_found"]).sum())
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    if n_no_metal:
        print(
            f"  [note] {n_no_metal}/{len(df)} structure(s) had no locatable metal point "
            "(no DDXXD / no coordinating O; e.g. OSC class-12 outlier) -> features NaN."
        )
    return df
