from __future__ import annotations

"""Active-site pocket descriptors for class-I TPS designs (fold-agnostic).

This metric characterises the *catalytic cavity* of a terpene-synthase design with
two independent cavity-detection engines, both anchored on the same physical point
— the active-site **metal point**: the centroid of the DDXXD + NSE/DTE coordinating
side-chain oxygens (exactly the centroid that ``active_site_geometry`` computes;
that module is reused here as the single source of truth for the motif logic).

Why anchor on the metal point? A structure can have many surface pockets; the one
that matters for TPS catalysis is the cavity that holds the trinuclear Mg2+/Mn2+
cluster and the prenyl-diphosphate substrate. By selecting the detected pocket that
encloses / is nearest the metal point we report descriptors for the *catalytic*
pocket specifically, not the largest or top-ranked one.

Two engines, run per structure:

* **fpocket** (geometric, Voronoi alpha-spheres). We pick the fpocket pocket whose
  alpha-spheres enclose the metal point (point inside the alpha-sphere cloud) or,
  failing that, the pocket whose alpha-sphere centroid is nearest the metal point.
  Reported columns:
    - ``catalytic_pocket_volume``      (A^3) — headline TPS signal; the molecular-
      ruler cavity volume tracks product chain length.
    - ``pocket_hydrophobicity``        — fpocket "Hydrophobicity Score".
    - ``pocket_enclosure``             — buriedness; fpocket "Polar Sasa /
      Apolar Sasa"... actually the dedicated descriptor is the proportion of
      apolar alpha spheres / mean local hydrophobic density; we surface fpocket's
      own enclosure-related score (see ``_FPOCKET_FIELDS``).
    - ``pocket_n_alpha_spheres``       — number of alpha spheres in the pocket.
    - ``pocket_total_sasa``            — total SASA of the pocket (A^2), if parsed.
    - ``pocket_depth``                 — pocket "max dist between apolar..." proxy;
      we surface fpocket's mean alpha-sphere radius as a depth proxy when no
      explicit depth field exists.
    - ``pocket_sasa_per_volume``       — DERIVED specific surface area =
      ``pocket_total_sasa / catalytic_pocket_volume`` (A^-1): a standard
      shape/compactness descriptor (higher = more surface per unit volume). NOTE it
      is size-dependent (~1/radius), so it co-varies with cavity size rather than
      isolating shape; a size-free shape factor (sphericity) would need a matched
      cavity surface+volume, which fpocket's SASA (lining-atom) and volume
      (alpha-sphere) are NOT. NaN when either input is missing / volume <= 0.
* **P2Rank** (machine-learned ligandability cross-check). We take the P2Rank
  predicted pocket nearest the metal point and report:
    - ``p2rank_catalytic_site_score``  — its ligandability ``score``.
    - ``p2rank_catalytic_pocket_rank`` — its 1-based rank among P2Rank's pockets
      (1 = P2Rank's top-ranked pocket). Independent confirmation that the catalytic
      cavity is a real ligand-binding site.

RAW numbers only — no bands (the reference-stats pipeline supplies the natural
band). When NO detected pocket coincides with the metal point (within a generous
cutoff), the engine's columns are NaN — itself a meaningful red flag and recorded
as such. When either motif is absent the metal point can't be located and ALL
descriptor columns are NaN.

CLI/dir conventions (af3-vs-flat detection, ID = filename stem,
``<structs_dir>_pocket_descriptors.csv`` naming) mirror ``plddt.py`` /
``active_site_geometry.py``.
"""

import argparse
import glob
import math
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBIO import PDBIO

# Reuse the sibling active-site-geometry machinery (motif localization + oxygen
# gathering) so the metal point is defined identically to that tool.
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from structure_metrics.active_site_geometry import (  # noqa: E402
    metal_point as _cage_metal_point,
    structure_sequence_residues_atoms,
)

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

COLUMNS = [
    "ID",
    "metal_point_found",
    # fpocket (geometric)
    "catalytic_pocket_volume",
    "pocket_hydrophobicity",
    "pocket_enclosure",
    "pocket_n_alpha_spheres",
    "pocket_total_sasa",
    "pocket_depth",
    "pocket_sasa_per_volume",
    "fpocket_catalytic_pocket_found",
    # P2Rank (ML cross-check)
    "p2rank_catalytic_site_score",
    "p2rank_catalytic_pocket_rank",
    "p2rank_catalytic_pocket_found",
    "n_residues",
]

# A detected pocket "coincides" with the metal point if its representative point is
# within this distance (A). Catalytic pockets are large; the alpha-sphere cloud /
# P2Rank centre routinely sits several A from the exact oxygen centroid, so use a
# generous cutoff. Beyond this we treat the catalytic pocket as not found (NaN).
METAL_POINT_CUTOFF_A = 12.0

# Mapping of fpocket per-pocket descriptor labels (as printed in fpocket 4.0's
# *_info.txt block) to our column names, matched case-insensitively as a substring
# of the label text BEFORE the ':'. Order matters: more specific needles first so
# e.g. "Volume score" doesn't shadow "Volume". The label text in fpocket 4.0 is e.g.
# "Volume :", "Hydrophobicity score:", "Total SASA :", "Number of Alpha Spheres :",
# "Mean local hydrophobic density :".
_FPOCKET_FIELDS = (
    ("number of alpha spheres", "pocket_n_alpha_spheres"),
    ("total sasa", "pocket_total_sasa"),
    ("hydrophobicity score", "pocket_hydrophobicity"),
    # enclosure / buriedness proxy: density of hydrophobic alpha-sphere neighbours.
    ("mean local hydrophobic density", "pocket_enclosure"),
    ("volume", "catalytic_pocket_volume"),
)


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def _nan_result(n_residues: int = 0) -> Dict[str, float]:
    return {
        "metal_point_found": False,
        "catalytic_pocket_volume": np.nan,
        "pocket_hydrophobicity": np.nan,
        "pocket_enclosure": np.nan,
        "pocket_n_alpha_spheres": np.nan,
        "pocket_total_sasa": np.nan,
        "pocket_depth": np.nan,
        "pocket_sasa_per_volume": np.nan,
        "fpocket_catalytic_pocket_found": False,
        "p2rank_catalytic_site_score": np.nan,
        "p2rank_catalytic_pocket_rank": np.nan,
        "p2rank_catalytic_pocket_found": False,
        "n_residues": n_residues,
    }


# --------------------------------------------------------------------------- #
# Metal point                                                                 #
# --------------------------------------------------------------------------- #
def metal_point(structure_path: str) -> Tuple[Optional[np.ndarray], int]:
    """Active-site metal point for one structure, via the CANONICAL relaxed definition
    in ``active_site_geometry.metal_point``: the centroid of the DDXXD (+ NSE/DTE when
    that motif is also matched) coordinating side-chain oxygens. Returns
    ``(point, n_residues)``; ``point`` is None when DDXXD is absent or no coordinating
    oxygen is found. For a both-motif structure this is identical to the prior
    both-motif centroid (the relaxation is a strict superset); it ADDS coverage for
    DDXXD-only real TPS (e.g. TEAS/5EAT) that the both-motif version left NaN."""
    sequence, residues, _ = structure_sequence_residues_atoms(structure_path)
    return _cage_metal_point(sequence, residues), len(sequence)


def _ensure_pdb(structure_path: str, workdir: str) -> str:
    """fpocket and P2Rank both want a PDB file. If the input is a .cif, convert it
    to PDB with Biopython into ``workdir`` and return that path; otherwise return a
    copy of the .pdb in workdir (so all tool output lands in the scratch dir)."""
    base = os.path.basename(structure_path)
    stem = os.path.splitext(base)[0]
    out_pdb = os.path.join(workdir, stem + ".pdb")
    if structure_path.lower().endswith((".cif", ".mmcif")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = _CIF_PARSER.get_structure("s", structure_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_pdb)
    else:
        shutil.copyfile(structure_path, out_pdb)
    return out_pdb


# --------------------------------------------------------------------------- #
# fpocket                                                                     #
# --------------------------------------------------------------------------- #
def _parse_fpocket_alpha_spheres(vert_pqr: str) -> np.ndarray:
    """Alpha-sphere centre coordinates (N,3) from an fpocket ``pocketN_vert.pqr``.
    fpocket writes each Voronoi vertex (alpha-sphere centre) as an ATOM record in
    residue 'STP', PDB-style coordinate columns 31-54. (The sibling
    ``pocketN_atm.pdb`` holds the pocket's PROTEIN atoms, not the alpha spheres.)"""
    pts: List[List[float]] = []
    try:
        with open(vert_pqr) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "STP":
                    pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    except (OSError, ValueError):
        pass
    return np.asarray(pts, dtype=float) if pts else np.empty((0, 3))


def _parse_fpocket_info(info_path: str) -> Dict[int, Dict[str, float]]:
    """Parse fpocket's ``*_info.txt`` into {pocket_number: {descriptor: value}}."""
    pockets: Dict[int, Dict[str, float]] = {}
    current: Optional[int] = None
    try:
        with open(info_path) as fh:
            for raw in fh:
                line = raw.strip()
                if line.lower().startswith("pocket"):
                    # "Pocket 1 :"
                    try:
                        current = int(line.split()[1])
                        pockets[current] = {}
                    except (IndexError, ValueError):
                        current = None
                    continue
                if current is None or ":" not in line:
                    continue
                label, _, val = line.partition(":")
                label = label.strip().lower()
                try:
                    value = float(val.strip())
                except ValueError:
                    continue
                # "volume" must be the exact descriptor, not "volume score".
                if label == "volume":
                    pockets[current]["catalytic_pocket_volume"] = value
                    continue
                for needle, col in _FPOCKET_FIELDS:
                    if col == "catalytic_pocket_volume":
                        continue  # handled above (exact match) to avoid "volume score"
                    if needle in label:
                        pockets[current][col] = value
                        break
    except OSError:
        pass
    return pockets


def run_fpocket(pdb_path: str, workdir: str, fpocket_bin: str = "fpocket") -> Optional[str]:
    """Run fpocket on ``pdb_path``; return the path to its ``*_out`` output dir, or
    None on failure. fpocket writes ``<stem>_out/`` next to the input PDB."""
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    try:
        subprocess.run(
            [fpocket_bin, "-f", pdb_path],
            cwd=workdir,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"  [warn] fpocket failed on {stem}: {exc}")
        return None
    out_dir = os.path.join(workdir, stem + "_out")
    return out_dir if os.path.isdir(out_dir) else None


def _point_inside_cloud(point: np.ndarray, cloud: np.ndarray, pad: float = 3.0) -> bool:
    """Crude enclosure test: is ``point`` inside the axis-aligned bounding box of the
    alpha-sphere cloud, padded by ``pad`` (A)? Alpha spheres tile the cavity, so a
    point within their padded extent is effectively enclosed by the pocket."""
    if len(cloud) == 0:
        return False
    lo = cloud.min(axis=0) - pad
    hi = cloud.max(axis=0) + pad
    return bool(np.all(point >= lo) and np.all(point <= hi))


def fpocket_catalytic(
    pdb_path: str, point: np.ndarray, workdir: str, fpocket_bin: str = "fpocket"
) -> Dict[str, float]:
    """Run fpocket and return the catalytic pocket's descriptors (the pocket whose
    alpha spheres enclose / are nearest the metal point). All-NaN/False when fpocket
    finds no pocket near the metal point."""
    res: Dict[str, float] = {
        "catalytic_pocket_volume": np.nan,
        "pocket_hydrophobicity": np.nan,
        "pocket_enclosure": np.nan,
        "pocket_n_alpha_spheres": np.nan,
        "pocket_total_sasa": np.nan,
        "pocket_depth": np.nan,
        "fpocket_catalytic_pocket_found": False,
    }
    out_dir = run_fpocket(pdb_path, workdir, fpocket_bin=fpocket_bin)
    if out_dir is None:
        return res

    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    info = _parse_fpocket_info(os.path.join(out_dir, stem + "_info.txt"))
    pockets_dir = os.path.join(out_dir, "pockets")
    pocket_verts = sorted(glob.glob(os.path.join(pockets_dir, "pocket*_vert.pqr")))

    best_num: Optional[int] = None
    best_dist = math.inf
    best_enclosed = False
    best_radius = np.nan
    for vpqr in pocket_verts:
        name = os.path.basename(vpqr)  # pocket12_vert.pqr
        try:
            num = int(name.replace("pocket", "").split("_")[0])
        except ValueError:
            continue
        cloud = _parse_fpocket_alpha_spheres(vpqr)
        if len(cloud) == 0:
            continue
        dists = np.sqrt(((cloud - point) ** 2).sum(axis=1))
        dmin = float(dists.min())
        enclosed = _point_inside_cloud(point, cloud)
        # Prefer an enclosing pocket; among those (or among non-enclosing if none
        # enclose) pick the nearest alpha sphere to the metal point.
        better = (enclosed and not best_enclosed) or (
            enclosed == best_enclosed and dmin < best_dist
        )
        if better:
            best_enclosed = enclosed
            best_dist = dmin
            best_num = num
            # mean alpha-sphere radius as a depth proxy (cavity scale)
            best_radius = float(np.mean(np.sqrt(((cloud - cloud.mean(axis=0)) ** 2).sum(axis=1))))

    if best_num is None or best_dist > METAL_POINT_CUTOFF_A:
        return res

    res["fpocket_catalytic_pocket_found"] = True
    desc = info.get(best_num, {})
    for col in ("catalytic_pocket_volume", "pocket_hydrophobicity",
                "pocket_enclosure", "pocket_n_alpha_spheres", "pocket_total_sasa"):
        if col in desc:
            res[col] = desc[col]
    res["pocket_depth"] = best_radius
    return res


# --------------------------------------------------------------------------- #
# P2Rank                                                                      #
# --------------------------------------------------------------------------- #
def run_p2rank(pdb_path: str, workdir: str, p2rank_bin: str) -> Optional[str]:
    """Run P2Rank ``predict -f <pdb>`` into ``workdir/p2rank_out``; return the path
    to the predictions CSV ``<name>.pdb_predictions.csv`` or None on failure."""
    out_dir = os.path.join(workdir, "p2rank_out")
    try:
        subprocess.run(
            [p2rank_bin, "predict", "-f", pdb_path, "-o", out_dir],
            cwd=workdir,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        stderr = getattr(exc, "stderr", b"") or b""
        print(f"  [warn] P2Rank failed on {os.path.basename(pdb_path)}: {exc} "
              f"{stderr.decode(errors='replace')[:300]}")
        return None
    base = os.path.basename(pdb_path)
    pred = os.path.join(out_dir, base + "_predictions.csv")
    if os.path.isfile(pred):
        return pred
    # P2Rank may strip the extension in the output name; glob as a fallback.
    hits = glob.glob(os.path.join(out_dir, "*_predictions.csv"))
    return hits[0] if hits else None


def _parse_p2rank_predictions(pred_csv: str) -> List[Dict[str, float]]:
    """Parse P2Rank predictions CSV into a list of pockets (in file order, i.e. by
    descending score = rank). Columns are whitespace-padded; we read center_x/y/z
    and score."""
    try:
        df = pd.read_csv(pred_csv, skipinitialspace=True)
    except Exception:
        return []
    df.columns = [c.strip() for c in df.columns]
    cols = {c.lower(): c for c in df.columns}
    needed = ["center_x", "center_y", "center_z", "score"]
    if not all(n in cols for n in needed):
        return []
    pockets: List[Dict[str, float]] = []
    for i, row in df.iterrows():
        pockets.append(
            {
                "rank": int(i) + 1,  # file order = rank (P2Rank sorts by score desc)
                "center": np.array(
                    [float(row[cols["center_x"]]), float(row[cols["center_y"]]),
                     float(row[cols["center_z"]])]
                ),
                "score": float(row[cols["score"]]),
            }
        )
    return pockets


def p2rank_catalytic(
    pdb_path: str, point: np.ndarray, workdir: str, p2rank_bin: str
) -> Dict[str, float]:
    """Run P2Rank and return the predicted pocket nearest the metal point: its
    ligandability score and 1-based rank. NaN/False when P2Rank finds no pocket
    within the cutoff."""
    res: Dict[str, float] = {
        "p2rank_catalytic_site_score": np.nan,
        "p2rank_catalytic_pocket_rank": np.nan,
        "p2rank_catalytic_pocket_found": False,
    }
    pred_csv = run_p2rank(pdb_path, workdir, p2rank_bin)
    if pred_csv is None:
        return res
    pockets = _parse_p2rank_predictions(pred_csv)
    if not pockets:
        return res
    best = min(pockets, key=lambda p: float(np.sqrt(((p["center"] - point) ** 2).sum())))
    dist = float(np.sqrt(((best["center"] - point) ** 2).sum()))
    if dist > METAL_POINT_CUTOFF_A:
        return res
    res["p2rank_catalytic_site_score"] = best["score"]
    res["p2rank_catalytic_pocket_rank"] = int(best["rank"])
    res["p2rank_catalytic_pocket_found"] = True
    return res


# --------------------------------------------------------------------------- #
# Per-structure + dir driver                                                  #
# --------------------------------------------------------------------------- #
def pocket_descriptors(
    structure_path: str,
    *,
    fpocket_bin: str = "fpocket",
    p2rank_bin: Optional[str] = None,
) -> Dict[str, float]:
    """fpocket + (optional) P2Rank catalytic-pocket descriptors for one structure,
    anchored on the active-site metal point. P2Rank is skipped (its columns stay
    NaN) when ``p2rank_bin`` is None."""
    point, n_residues = metal_point(structure_path)
    result = _nan_result(n_residues)
    if point is None:
        return result  # metal point not locatable -> all descriptors NaN
    result["metal_point_found"] = True

    workdir = tempfile.mkdtemp(prefix="pocket_")
    try:
        pdb_path = _ensure_pdb(structure_path, workdir)
        result.update(fpocket_catalytic(pdb_path, point, workdir, fpocket_bin=fpocket_bin))
        if p2rank_bin:
            result.update(p2rank_catalytic(pdb_path, point, workdir, p2rank_bin))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    # Derived size-relative descriptor: specific surface area = Total SASA / Volume
    # (A^-1). Size-dependent (∝1/radius), but a standard pocket shape/compactness
    # descriptor. NaN when either fpocket quantity is missing or volume is non-positive.
    v = result.get("catalytic_pocket_volume", np.nan)
    a = result.get("pocket_total_sasa", np.nan)
    result["pocket_sasa_per_volume"] = (
        float(a / v) if (np.isfinite(v) and np.isfinite(a) and v > 0) else np.nan
    )
    return result


def _collect_structures(structs_dir: str) -> Tuple["OrderedDict[str, str]", str]:
    """Map ID -> structure file, auto-detecting layout (mirrors plddt.py /
    active_site_geometry.py): af3 ``af_output`` dir (``<job>/<job>_model.cif``;
    ID = job name) takes precedence; otherwise a flat dir of .pdb/.cif (ID =
    filename stem; .pdb wins on tie)."""
    af3: Dict[str, str] = {}
    try:
        entries = sorted(os.listdir(structs_dir))
    except FileNotFoundError:
        entries = []
    for entry in entries:
        sub = os.path.join(structs_dir, entry)
        model = os.path.join(sub, entry + "_model.cif")
        if os.path.isdir(sub) and os.path.isfile(model):
            af3[entry] = model
    if af3:
        return OrderedDict(sorted(af3.items())), "af3"

    chosen: Dict[str, str] = {}
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, f"*{ext}"))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return OrderedDict(sorted(chosen.items())), "flat"


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_pocket_descriptors.csv")


def pocket_descriptors_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    fpocket_bin: str = "fpocket",
    p2rank_bin: Optional[str] = None,
) -> pd.DataFrame:
    """fpocket (+ optional P2Rank) catalytic-pocket descriptors for every structure
    in ``structs_dir``; CSV keyed by ID. The catalytic pocket is the detected pocket
    enclosing / nearest the active-site metal point (DDXXD+NSE/DTE oxygen centroid).
    """
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")
    if p2rank_bin:
        print(f"Using P2Rank: {p2rank_bin}")
    else:
        print("[note] no P2Rank binary given (P2RANK_PATH unset); "
              "p2rank_* columns will be NaN.")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = pocket_descriptors(path, fpocket_bin=fpocket_bin, p2rank_bin=p2rank_bin)
        except Exception as exc:  # malformed/unparsable -> NaN row, keep going
            print(f"  [warn] failed on {os.path.basename(path)}: {exc}")
            stats = _nan_result(0)
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        print(f"  processed {i}/{n}: {stem}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)

    n_no_metal = int((~df["metal_point_found"]).sum())
    n_no_pocket = int((df["metal_point_found"] & ~df["fpocket_catalytic_pocket_found"]).sum())
    print(f"Wrote {len(df)} rows to {save_path}"
          + (f" ({n_failed} unparsable)" if n_failed else ""))
    if n_no_metal:
        print(f"  [flag] {n_no_metal}/{len(df)} structures: metal point not locatable "
              "(motif absent / no coordinating O) -> all descriptors NaN.")
    if n_no_pocket:
        print(f"  [flag] {n_no_pocket}/{len(df)} structures: metal point found but NO "
              "fpocket pocket coincides with it (red flag).")
    return df
