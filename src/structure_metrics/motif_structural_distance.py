from __future__ import annotations

"""3D (Angstrom) distance between the two class-I TPS metal-binding motifs.

Fold-agnostic: operates on a directory of structures (.pdb/.cif), so it works for
BOTH AlphaFold and ESMFold output with the same code. For each structure we

1. read the residue chain with Biopython (reusing plddt.py's parser selection and
   af_output/flat auto-detection),
2. derive the 1-letter sequence from the residues,
3. run the SHARED motif-localization helper on that sequence to find the
   DDXXD-family and NSE/DTE motifs, and
4. take the CA coordinates of the motif residues and compute the inter-motif 3D
   distance.

The two Mg/Mn ions sit in the cleft BETWEEN these motifs, so the distance between
the centroids of the two motifs' metal-coordinating CA atoms approximates the span
of the active-site metal cluster. We report:

* ``motif_centroid_distance`` — distance between the centroid of the DDXXD
  coordinating-residue CA atoms and the centroid of the NSE/DTE coordinating-
  residue CA atoms (the primary metric);
* ``motif_min_ca_distance`` — the minimum CA-CA distance between any DDXXD
  coordinating residue and any NSE/DTE coordinating residue (closest approach).

Output is a CSV keyed by ``ID`` (filename stem, matching the other structure
tools). Distances are NaN when a motif isn't found in the structure-derived
sequence. ``n_residues`` reports the modelled residue count for context.
"""

import argparse
import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Reuse the sequence-metrics shared motif localization (single source of truth).
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from sequence_metrics.motif_localization import (  # noqa: E402
    DDXXD_COORDINATING_OFFSETS,
    NSE_DTE_COORDINATING_OFFSETS,
    coordinating_indices,
    locate_ddxxd,
    locate_nse_dte,
)

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

COLUMNS = ["ID", "motif_centroid_distance", "motif_min_ca_distance", "n_residues"]


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def _three_to_one(resname: str) -> str:
    """3-letter residue name -> 1-letter code, 'X' for anything non-standard."""
    try:
        return index_to_one(three_to_index(resname))
    except KeyError:
        return "X"


def structure_sequence_and_ca(structure_path: str) -> Tuple[str, List[Optional[np.ndarray]]]:
    """Return (1-letter sequence, list of CA xyz arrays) for the protein residues
    of the first model, in chain order. The two lists are index-aligned: position i
    in the sequence has CA coordinate ``ca_coords[i]`` (None if that residue lacks a
    CA). HETATM ligands/ions/water are skipped — matching plddt.residue_plddts."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))  # first model only (predicted structs write one)
    seq_chars: List[str] = []
    ca_coords: List[Optional[np.ndarray]] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM (ions/ligands/water)
                continue
            seq_chars.append(_three_to_one(residue.get_resname()))
            ca_coords.append(
                np.asarray(residue["CA"].get_coord(), dtype=float) if "CA" in residue else None
            )
    return "".join(seq_chars), ca_coords


def _coord_matrix(indices: List[int], ca_coords: List[Optional[np.ndarray]]) -> Optional[np.ndarray]:
    """Stack the CA coords for the given residue indices, dropping any missing
    CA. Returns an (n, 3) array, or None if no usable coordinate remains."""
    pts = [ca_coords[i] for i in indices if 0 <= i < len(ca_coords) and ca_coords[i] is not None]
    return np.vstack(pts) if pts else None


def motif_distances(structure_path: str) -> Dict[str, float]:
    """Centroid and min CA-CA distance between the DDXXD and NSE/DTE motifs of one
    structure; NaN distances when a motif is absent or has no usable CA."""
    sequence, ca_coords = structure_sequence_and_ca(structure_path)
    n_residues = len(sequence)
    result = {
        "motif_centroid_distance": np.nan,
        "motif_min_ca_distance": np.nan,
        "n_residues": n_residues,
    }

    ddxxd = locate_ddxxd(sequence)
    nse = locate_nse_dte(sequence)
    if ddxxd is None or nse is None:
        return result

    ddxxd_pts = _coord_matrix(coordinating_indices(ddxxd, DDXXD_COORDINATING_OFFSETS), ca_coords)
    nse_pts = _coord_matrix(coordinating_indices(nse, NSE_DTE_COORDINATING_OFFSETS), ca_coords)
    if ddxxd_pts is None or nse_pts is None:
        return result

    centroid_dist = float(np.linalg.norm(ddxxd_pts.mean(axis=0) - nse_pts.mean(axis=0)))
    # Pairwise CA-CA distances between the two motifs' coordinating residues.
    diffs = ddxxd_pts[:, None, :] - nse_pts[None, :, :]
    min_dist = float(np.sqrt((diffs ** 2).sum(axis=2)).min())

    result["motif_centroid_distance"] = centroid_dist
    result["motif_min_ca_distance"] = min_dist
    return result


def _collect_structures(structs_dir: str) -> Tuple["OrderedDict[str, str]", str]:
    """Map ID -> structure file, auto-detecting layout. Mirrors plddt._collect_structures:
    an AF3 ``af_output`` dir (per-job ``<job>/<job>_model.cif``; ID = job name) takes
    precedence; otherwise a flat dir of .pdb/.cif (ID = filename stem; .pdb wins on tie)."""
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
    import glob
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, f"*{ext}"))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return OrderedDict(sorted(chosen.items())), "flat"


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_motif_structural_distance.csv")


def motif_structural_distance_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """Inter-motif 3D distance for every structure in ``structs_dir``; CSV keyed by ID."""
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
            stats = motif_distances(path)
        except Exception as exc:  # malformed/unparsable -> NaN row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {"motif_centroid_distance": np.nan, "motif_min_ca_distance": np.nan,
                     "n_residues": 0}
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
