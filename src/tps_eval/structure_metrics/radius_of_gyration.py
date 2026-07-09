from __future__ import annotations

"""Radius-of-gyration / compactness structural metrics for TPS designs.

Fold-agnostic global shape descriptors computed from the Cα coordinates of a
single design chain (AlphaFold/ESMFold .pdb/.cif). These are RAW geometric
numbers only — no expected-Rg band, no compactness ratio, no reference
comparison (that comparison is done downstream by a separate reference-statistics
pipeline). Every structure yields one row; an unparsable structure yields a NaN
row and the run continues.

Metrics (columns, keyed by ``ID``):

* ``radius_of_gyration`` — Rg = sqrt(mean_i ||r_i - r_com||^2) over the Cα atoms,
  in Å. UNWEIGHTED (each Cα counts equally; not mass-weighted — Cα masses are
  identical anyway, and we deliberately don't pull in heavy/side-chain atoms, so
  the unweighted Cα Rg is the natural choice). r_com is the Cα centroid.
* ``asphericity`` — from the 3×3 gyration tensor of the Cα coordinates with sorted
  eigenvalues λ1≥λ2≥λ3: b = λ1 − (λ2+λ3)/2. 0 = spherical; large = elongated /
  prolate (e.g. a two-domain or extended fold). Same length unit² convention as
  the gyration-tensor eigenvalues (Å²).
* ``acylindricity`` — c = λ2 − λ3 (0 = cylindrically symmetric about the major
  axis). Cheap byproduct of the same eigendecomposition.
* ``principal_radius_1/2/3`` — sqrt(λ1), sqrt(λ2), sqrt(λ3) in Å, the principal
  radii of gyration along the tensor's principal axes (note
  Rg = sqrt(λ1+λ2+λ3)).
* ``n_residues`` — number of Cα atoms used.

CHAIN CHOICE: a TPS design is a single chain, but to be robust we use ALL protein
Cα atoms across every chain of the first model (HETATM ions/ligands/water are
skipped via the residue hetflag, exactly as in ``plddt.py``). For a single-chain
design this is simply that chain's Cα atoms; for an (unexpected) multi-chain model
the Rg is computed over the union of protein Cα — documented here so the number is
interpretable.
"""

import glob
import os
import warnings
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

# Output column order (ID first, then the raw geometric metrics).
COLUMNS = [
    "ID",
    "radius_of_gyration",
    "asphericity",
    "acylindricity",
    "principal_radius_1",
    "principal_radius_2",
    "principal_radius_3",
    "n_residues",
]


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def ca_coordinates(structure_path: str) -> np.ndarray:
    """(N, 3) array of Cα coordinates: the CA atom of every protein residue that
    has one, across all chains of the first model. HETATM ions/ligands/water are
    skipped (hetflag != ' '), mirroring ``plddt.py``. Works for .pdb and .cif."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))  # first model only (predicted structs write one)
    coords: List[np.ndarray] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM (ions/ligands/water)
                continue
            if "CA" in residue:
                coords.append(np.asarray(residue["CA"].get_coord(), dtype=float))
    return np.vstack(coords) if coords else np.empty((0, 3))


def gyration_metrics(coords: np.ndarray) -> Dict[str, float]:
    """Raw radius-of-gyration / shape metrics from a set of (N, 3) coordinates.

    Rg uses the unweighted definition sqrt(mean ||r_i - r_com||^2). The gyration
    tensor S = (1/N) Σ (r_i - r_com)(r_i - r_com)^T has eigenvalues λ1≥λ2≥λ3 with
    Rg^2 = λ1+λ2+λ3; asphericity = λ1 − (λ2+λ3)/2, acylindricity = λ2 − λ3, and
    the principal radii are sqrt(λk). NaN-everything for fewer than 1 Cα."""
    nan_result: Dict[str, float] = {
        "radius_of_gyration": np.nan,
        "asphericity": np.nan,
        "acylindricity": np.nan,
        "principal_radius_1": np.nan,
        "principal_radius_2": np.nan,
        "principal_radius_3": np.nan,
        "n_residues": 0,
    }
    if coords.shape[0] == 0:
        return nan_result

    com = coords.mean(axis=0)
    centered = coords - com
    n = centered.shape[0]

    rg = float(np.sqrt((centered ** 2).sum(axis=1).mean()))

    # Gyration tensor (1/N Σ outer products). Symmetric & PSD -> eigvalsh.
    tensor = (centered.T @ centered) / n
    eig = np.linalg.eigvalsh(tensor)  # ascending
    lam3, lam2, lam1 = float(eig[0]), float(eig[1]), float(eig[2])  # λ1≥λ2≥λ3

    return {
        "radius_of_gyration": rg,
        "asphericity": lam1 - 0.5 * (lam2 + lam3),
        "acylindricity": lam2 - lam3,
        "principal_radius_1": float(np.sqrt(max(lam1, 0.0))),
        "principal_radius_2": float(np.sqrt(max(lam2, 0.0))),
        "principal_radius_3": float(np.sqrt(max(lam3, 0.0))),
        "n_residues": int(n),
    }


def radius_of_gyration(structure_path: str) -> Dict[str, float]:
    """Compute the raw gyration/compactness metrics for one structure."""
    return gyration_metrics(ca_coordinates(structure_path))


def _collect_structures(structs_dir: str) -> Tuple["OrderedDict[str, str]", str]:
    """Map ID -> structure file, auto-detecting layout. Mirrors
    plddt/active_site_geometry: an AF3 ``af_output`` dir (per-job
    ``<job>/<job>_model.cif``; ID = job name) takes precedence; otherwise a flat
    dir of .pdb/.cif (ID = filename stem; .pdb wins on tie)."""
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
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_radius_of_gyration.csv")


def radius_of_gyration_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """Raw radius-of-gyration / compactness metrics for every structure in
    ``structs_dir``; CSV keyed by ID. Unparsable structures yield a NaN row."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = radius_of_gyration(path)
        except Exception as exc:  # malformed/unparsable -> NaN row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = gyration_metrics(np.empty((0, 3)))
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    return df
