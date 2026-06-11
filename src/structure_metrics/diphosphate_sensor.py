from __future__ import annotations

"""Diphosphate-sensor active-site metric for class-I TPS designs.

This rounds out the catalytic checklist beyond the Mg-binding carboxylates (the
``active_site_geometry`` carboxylate cage) and the aromatic pocket lining
(``aromatic_lining``). After the trinuclear Mg2+ cluster forms (chelated by the
DDXXD + NSE/DTE carboxylates), class-I terpene synthases use a set of *basic*
residues — a pyrophosphate-**sensor Arg**, often paired with a Tyr (the conserved
"RY pair"), plus additional Arg/Lys — to bind the substrate's diphosphate (PPi)
and help trigger ionization (the induced-fit "effector triad"). This tool checks
that a design actually PRESENTS those residues, positioned at the metal /
diphosphate locus.

The diphosphate sits between the metal cluster and these basic residues, so we use
the carboxylate-cage **metal point** as the geometric anchor and look for Arg/Lys
side-chain terminal nitrogens that lie within a cutoff of it AND point toward it.

For each structure we

1. parse it with Biopython (reusing the af_output/flat auto-detection + ID-stem
   convention via ``sdr_divergence.collect_structures`` and the per-residue
   ``ResidueInfo``),
2. locate the active-site **metal point** by REUSING the canonical
   ``active_site_geometry.metal_point`` — the centroid of the DDXXD (+ NSE/DTE when
   that motif is also matched) coordinating side-chain oxygens. DDXXD is required and
   NSE/DTE is an optional refinement: the shared NSE/DTE regex misses real TEAS (5EAT),
   so requiring BOTH motifs would make real references all-NaN. The relaxed set is a
   strict superset — identical to the both-motif centroid when both motifs match.
3. count **basic residues near the diphosphate site**: Arg/Lys whose side-chain
   terminal N atoms (Arg NE/NH1/NH2, Lys NZ) lie within ``--cutoff`` A (default 12)
   of the metal point AND point toward it (the residue's terminal-N is closer to the
   metal point than its Cα is — i.e. the basic head reaches in, rather than the side
   chain pointing away). Report ``n_diphosphate_basic_residues`` (+ ``n_arg`` /
   ``n_lys``).
4. detect the **RY pair**: an Arg counted in (3) that has a partner Tyr either
   adjacent in sequence (within +-2 residues) OR spatially close (Tyr OH within
   ``--ry_dist`` A, default 6, of that Arg's guanidinium centroid), the Tyr itself
   being near the site (OH within ``--cutoff`` of the metal point). Report
   ``has_RY_pair`` (bool) and ``n_RY_pairs``.

Columns (keyed by ``ID``):

* ``metal_point_found``      — bool; True when the carboxylate-cage anchor was located.
* ``n_diphosphate_basic_residues`` — Arg+Lys near & pointing at the site.
* ``n_arg`` / ``n_lys``      — the breakdown.
* ``has_RY_pair`` (bool) / ``n_RY_pairs`` — conserved sensor-Arg / Tyr pairs.
* ``n_residues``             — modelled residue count, for context.

NaN/0 when the site can't be located: ``metal_point_found`` is then False, the
counts are 0 and the flags False (a real flag, not a silent zero). Every structure
gets exactly one row.

APO CAVEAT: residue PRESENCE and COUNT are robust on apo open-state models; the
exact rotamer orientation is softer (the basic side chains swing in on substrate /
metal binding). We therefore treat presence + a COARSE direction check (terminal-N
nearer the metal point than the Cα) as the primary signal, not a tight geometric
gate. The spatial RY-pair criterion is likewise generous (6 A) for the same reason.
"""

import argparse
import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Reuse the sibling structure-metrics + specificity modules (single source of truth
# for structure loading, per-residue access, and the metal-point computation WITH
# the DDXXD-only NSE/DTE fallback). Do NOT re-encode any of these.
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from specificity.sdr_divergence import (  # noqa: E402
    ResidueInfo,
    collect_structures,
)
from structure_metrics.active_site_geometry import metal_point  # noqa: E402

# --------------------------------------------------------------------------- #
# Defaults (documented; apo-robust where possible)
# --------------------------------------------------------------------------- #
# Distance (A) from the metal point within which a basic residue's terminal N must
# lie to count as a diphosphate-site residue. The PPi sits between the metals and
# these residues; ~12 A comfortably spans the metal point -> PPi -> guanidinium/
# amine reach without pulling in bulk-surface basics.
DEFAULT_CUTOFF = 12.0

# Distance (A) for the SPATIAL RY-pair criterion: a Tyr OH within this of a sensor
# Arg's guanidinium centroid. Generous (apo rotamers are soft); the sequence-
# adjacency criterion (+-RY_SEQ_WINDOW) is the orientation-independent complement.
DEFAULT_RY_DIST = 6.0

# Sequence-adjacency window (residues) for the "RY pair adjacent in sequence" check.
RY_SEQ_WINDOW = 2

# Arg side-chain terminal nitrogens (the guanidinium) and Lys terminal nitrogen.
ARG_TERMINAL_N = ("NE", "NH1", "NH2")
LYS_TERMINAL_N = ("NZ",)
ARG_GUANIDINIUM = ("NE", "NH1", "NH2", "CZ")  # for the guanidinium centroid

COLUMNS = [
    "ID",
    "metal_point_found",
    "n_diphosphate_basic_residues",
    "n_arg",
    "n_lys",
    "has_RY_pair",
    "n_RY_pairs",
    "n_residues",
]


def _atom_coords(residue, atom_names) -> List[np.ndarray]:
    """Coords of the named atoms present in a Biopython residue."""
    out: List[np.ndarray] = []
    for name in atom_names:
        if name in residue:
            out.append(np.asarray(residue[name].get_coord(), dtype=float))
    return out


def _ca_coord(residue) -> Optional[np.ndarray]:
    if "CA" in residue:
        return np.asarray(residue["CA"].get_coord(), dtype=float)
    return None


def diphosphate_sensor_one(
    info: ResidueInfo,
    *,
    cutoff: float = DEFAULT_CUTOFF,
    ry_dist: float = DEFAULT_RY_DIST,
) -> Dict[str, object]:
    """Diphosphate-sensor metrics for one parsed structure.

    Locates the carboxylate-cage metal point (canonical relaxed
    ``active_site_geometry.metal_point``), counts the Arg/Lys whose terminal N atoms are
    within ``cutoff`` of it AND point toward it, and detects RY (sensor-Arg / Tyr)
    pairs. Returns the row dict (without ID)."""
    row: Dict[str, object] = {
        "metal_point_found": False,
        "n_diphosphate_basic_residues": 0,
        "n_arg": 0,
        "n_lys": 0,
        "has_RY_pair": False,
        "n_RY_pairs": 0,
        "n_residues": len(info.seq),
    }

    mp = metal_point(info.seq, info.residues)
    if mp is None:
        return row  # site not locatable -> all-zero / False, metal_point_found stays False
    row["metal_point_found"] = True

    cutoff2 = cutoff * cutoff

    # Pass 1: collect the basic residues near & pointing at the site. Record each
    # qualifying Arg's guanidinium centroid + its sequence index so the RY-pair pass
    # can pair it with a nearby Tyr.
    arg_hits: List[Dict[str, object]] = []  # {idx, guanidinium_centroid}
    n_arg = 0
    n_lys = 0
    for i, residue in enumerate(info.residues):
        resname = residue.get_resname()
        if resname == "ARG":
            terminal = _atom_coords(residue, ARG_TERMINAL_N)
        elif resname == "LYS":
            terminal = _atom_coords(residue, LYS_TERMINAL_N)
        else:
            continue
        if not terminal:
            continue
        terminal_arr = np.vstack(terminal)
        d2 = ((terminal_arr - mp) ** 2).sum(axis=1)
        nearest2 = float(d2.min())
        if nearest2 > cutoff2:
            continue
        # Coarse direction check: the basic head must reach TOWARD the metal point —
        # its nearest terminal N is closer to the metal point than its Cα is. Soft on
        # apo open-state models, hence "coarse" (see APO CAVEAT in the module docstring).
        ca = _ca_coord(residue)
        if ca is not None:
            ca_d2 = float(((ca - mp) ** 2).sum())
            if nearest2 > ca_d2:
                continue  # side chain points away from the site
        if resname == "ARG":
            n_arg += 1
            guan = _atom_coords(residue, ARG_GUANIDINIUM)
            centroid = np.vstack(guan).mean(axis=0) if guan else terminal_arr.mean(axis=0)
            arg_hits.append({"idx": i, "centroid": centroid})
        else:
            n_lys += 1

    row["n_arg"] = n_arg
    row["n_lys"] = n_lys
    row["n_diphosphate_basic_residues"] = n_arg + n_lys

    # Pass 2: RY pairs. For each qualifying sensor Arg, look for a partner Tyr that is
    # (a) adjacent in sequence (within +-RY_SEQ_WINDOW), OR (b) spatially close (its OH
    # within ry_dist of the Arg's guanidinium centroid) AND itself near the site (OH
    # within cutoff of the metal point). Each Arg pairs at most once (counted once).
    n_ry = 0
    for hit in arg_hits:
        ai = int(hit["idx"])
        centroid = hit["centroid"]
        paired = False
        for j, residue in enumerate(info.residues):
            if residue.get_resname() != "TYR":
                continue
            oh = _atom_coords(residue, ("OH",))
            adjacent = abs(j - ai) <= RY_SEQ_WINDOW
            spatial = False
            if oh:
                oh_xyz = oh[0]
                near_arg = float(((oh_xyz - centroid) ** 2).sum()) <= ry_dist * ry_dist
                near_site = float(((oh_xyz - mp) ** 2).sum()) <= cutoff2
                spatial = near_arg and near_site
            if adjacent or spatial:
                paired = True
                break
        if paired:
            n_ry += 1

    row["n_RY_pairs"] = n_ry
    row["has_RY_pair"] = bool(n_ry > 0)
    return row


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_diphosphate_sensor.csv")


def diphosphate_sensor_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    cutoff: float = DEFAULT_CUTOFF,
    ry_dist: float = DEFAULT_RY_DIST,
) -> pd.DataFrame:
    """Diphosphate-sensor metric for every structure in ``structs_dir``; CSV keyed by ID.

    Mirrors the other structure-branch tools: auto-detects an AF3 ``af_output`` dir
    vs a flat dir of .pdb/.cif (via ``sdr_divergence.collect_structures``), writes
    ``<structs_dir>_diphosphate_sensor.csv`` by default, one row per structure."""
    structures, mode = collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    rows: List[Dict[str, object]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                info = ResidueInfo(path)
            stats = diphosphate_sensor_one(info, cutoff=cutoff, ry_dist=ry_dist)
        except Exception as exc:  # malformed/unparsable -> NaN/0 row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {
                "metal_point_found": False,
                "n_diphosphate_basic_residues": 0,
                "n_arg": 0,
                "n_lys": 0,
                "has_RY_pair": False,
                "n_RY_pairs": 0,
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
    n_found = int(df["metal_point_found"].sum())
    print(
        f"Wrote {len(df)} rows to {save_path} "
        f"({n_found}/{len(df)} with a locatable metal point"
        + (f", {n_failed} unparsable" if n_failed else "")
        + ")."
    )
    return df
