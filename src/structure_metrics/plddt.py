from __future__ import annotations

import glob
import os
import warnings
from collections import OrderedDict
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# AlphaFold stores per-residue pLDDT (0-100) in the B-factor field of both PDB
# (cols 61-66) and mmCIF (_atom_site.B_iso_or_equiv). pLDDT >= 70 is the common
# "confident" cutoff; >= 90 is "very high".
#
# IMPORTANT: this reads the B-factor field as pLDDT, which is only valid for
# predicted structures (AlphaFold / ESMFold / ...). For EXPERIMENTAL structures
# (PDB depositions, X-ray/cryo-EM) the B-factor is the crystallographic
# temperature factor, NOT a confidence score — the numbers it produces for those
# are meaningless. Point this at a directory of predicted structures only.
CONFIDENT_THRESHOLD = 70.0

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

# Output column order (ID first, then the filtration metrics).
COLUMNS = ["ID", "mean_plddt", "median_plddt", "min_plddt", "frac_plddt_confident", "n_residues"]


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def residue_plddts(structure_path: str) -> List[float]:
    """Per-residue pLDDT for a structure: the CA B-factor of every residue that
    has one, across all chains of the first model. Works for both .pdb and .cif."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))  # first model only (AlphaFold writes one)
    plddts: List[float] = []
    for chain in model:
        for residue in chain:
            # Standard polymer (protein/nucleic) residues only. Skip HETATM
            # ligands/ions/water (hetflag != " "): otherwise a kept ion whose atom
            # is named "CA" (e.g. a calcium ion) would be miscounted as a residue's
            # pLDDT. (No-op for protein-only structures.)
            if residue.id[0] != " ":
                continue
            if "CA" in residue:
                plddts.append(float(residue["CA"].get_bfactor()))
    return plddts


def summarize(plddts: List[float], confident_threshold: float = CONFIDENT_THRESHOLD) -> Dict[str, float]:
    """Per-structure summary stats used as filtration criteria."""
    arr = np.asarray(plddts, dtype=float)
    if arr.size == 0:
        return {
            "mean_plddt": np.nan,
            "median_plddt": np.nan,
            "min_plddt": np.nan,
            "frac_plddt_confident": np.nan,
            "n_residues": 0,
        }
    return {
        "mean_plddt": float(arr.mean()),
        "median_plddt": float(np.median(arr)),
        "min_plddt": float(arr.min()),
        "frac_plddt_confident": float((arr >= confident_threshold).mean()),
        "n_residues": int(arr.size),
    }


def _collect_structures(structs_dir: str):
    """Map ID -> structure file, auto-detecting the input layout. Returns
    (OrderedDict[id -> path], mode).

    Two layouts are supported:

    * "af3"  — an AlphaFold3 ``af_output`` directory: one subfolder per job,
      each containing the top-ranked ``<job>/<job>_model.cif`` (whose B-factor
      column IS pLDDT). This is the AUTHORITATIVE source. ID = job subfolder
      name (== the structs/ stem in the canonical pipeline). Note: structs/*.pdb
      extracted by the OLD Open Babel cif->pdb converter have their B-factor
      zeroed; vendor/cif_to_pdb was patched (Biopython) to preserve pLDDT, so
      freshly-extracted structs/ now carry it — but re-extract older runs or use
      af_output if you see all-zero pLDDT.

    * "flat" — a directory of structure files (``.pdb``/``.cif``) whose B-factor
      already holds pLDDT (e.g. AF2/ColabFold output, or any structures known to
      carry per-residue confidence). If both .pdb and .cif exist for an ID the
      .pdb wins (same values, simpler parser). ID = filename stem.
    """
    # AF3 af_output layout takes precedence when detected.
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

    # Flat directory of structure files. Order matters: later wins, so .pdb last.
    chosen: Dict[str, str] = {}
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, f"*{ext}"))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return OrderedDict(sorted(chosen.items())), "flat"


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_plddt.csv")


def extract_plddt_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    confident_threshold: float = CONFIDENT_THRESHOLD,
) -> pd.DataFrame:
    """Extract per-structure pLDDT summaries for every structure in `structs_dir`
    and write a CSV (keyed by ID) usable as a filtration criterion alongside the
    other tps_eval metrics."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} "
            "(expected an AlphaFold3 af_output dir with <job>/<job>_model.cif "
            "subfolders, or a flat dir of .pdb/.cif files)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = summarize(residue_plddts(path), confident_threshold)
        except Exception as exc:  # malformed/unparsable structure -> NaN row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = summarize([], confident_threshold)
            n_failed += 1
        stats["ID"] = stem
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    # All-zero pLDDT on residues that DID parse means the B-factor was never
    # populated — typically structs/*.pdb extracted by the old Open Babel
    # cif->pdb converter. Flag it loudly so the result isn't trusted as real pLDDT.
    n_zero = int(((df["n_residues"] > 0) & (df["mean_plddt"] == 0.0)).sum())
    if n_zero:
        print(
            f"  [warn] {n_zero}/{len(df)} structures have all-zero pLDDT — their B-factor "
            "looks unpopulated (e.g. structs/ extracted by the old obabel cif->pdb "
            "converter). Re-extract with the patched vendor/cif_to_pdb, or run on the "
            "AF3 af_output directory instead."
        )

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    return df
