"""Global fold confidence (pTM / ipTM) from the fold-time ``<ID>_pae.npz`` files.

pTM (predicted TM-score) is AlphaFold/ESMFold's WHOLE-FOLD confidence scalar: a
single number in [0, 1] (higher = better) estimating the TM-score between the
prediction and the (unknown) true structure. Unlike pLDDT — which is PER-RESIDUE
local confidence — pTM judges the GLOBAL fold / overall topology, so it complements
pLDDT as a filtration criterion (a design can be locally confident everywhere yet
have a globally uncertain fold). Well-folded TPS land around pTM ~0.7-0.9.

ipTM (interface pTM) is the multi-chain analogue (confidence in inter-chain
interfaces). TPS designs here are single-chain, so ipTM is usually absent (NaN); it
is emitted as a column only when at least one npz carries a non-NaN ipTM.

This tool does NOT re-fold and does NOT parse structures: the scalars are saved at
fold time into the SAME ``<ID>_pae.npz`` that interdomain_pae consumes (written by
esmfold.py / scripts/alphafold/extract_pae.py). So it just reads ``--pae_dir`` and
reduces to one row per ID. Output is the dir-keyed CSV ``<structs_dir>_global_confidence.csv``
(mirroring plddt's naming) keyed by ``ID`` (= the npz filename stem, == structure
stem), with columns ``ptm`` (always) and ``iptm`` (only if present anywhere). NaN
when the npz is missing or carries no ptm. Raw numbers only — "natural TPS" bands
are computed separately by the reference-stats pipeline.

The PAE npz schema (one schema, both folders — see esmfold.py / extract_pae.py):
  * ``pae``         : float32 (L, L), Angstrom.
  * ``residue_ids`` : int32 (L,), the PDB author residue number of each PAE row/col.
  * ``n_residues``  : int scalar (== L).
  * ``source``      : str, 'esmfold' | 'alphafold3'.
  * ``ptm``         : float32 scalar, global fold confidence (0-1; NaN if absent).
  * ``iptm``        : float32 scalar, interface pTM (multi-chain only; NaN if absent).
"""

from __future__ import annotations

import glob
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Base (always present) output columns; iptm is appended only when present anywhere.
BASE_COLUMNS = ["ID", "ptm"]
IPTM_COLUMN = "iptm"


def _pae_path(pae_dir: str, design_id: str) -> str:
    return os.path.join(pae_dir, f"{design_id}_pae.npz")


def _scalar_or_nan(npz, key: str) -> float:
    """Read a 0-d/1-element scalar field from an npz, NaN if absent/unreadable."""
    if key not in npz.files:
        return float("nan")
    try:
        return float(np.asarray(npz[key]).reshape(-1)[0])
    except Exception:  # noqa: BLE001 - odd shape/dtype -> treat as missing
        return float("nan")


def load_global_confidence(pae_path: str) -> Dict[str, float]:
    """Read ``ptm`` and ``iptm`` from one ``<ID>_pae.npz``.

    Returns ``{"ptm": float, "iptm": float}`` (NaN for a missing file or a missing
    field). Reading the scalars does NOT load the (L, L) PAE matrix into the result.
    """
    if not os.path.isfile(pae_path):
        return {"ptm": float("nan"), IPTM_COLUMN: float("nan")}
    with np.load(pae_path, allow_pickle=False) as npz:
        return {
            "ptm": _scalar_or_nan(npz, "ptm"),
            IPTM_COLUMN: _scalar_or_nan(npz, IPTM_COLUMN),
        }


def _ids_from_pae_dir(pae_dir: str) -> List[str]:
    """Every ``<ID>_pae.npz`` stem in `pae_dir` (ID = filename minus ``_pae.npz``)."""
    ids = []
    for path in sorted(glob.glob(os.path.join(pae_dir, "*_pae.npz"))):
        ids.append(os.path.basename(path)[: -len("_pae.npz")])
    return ids


def _default_save_path(structs_dir: str) -> str:
    """Dir-keyed CSV sibling of the structs dir, mirroring plddt's naming."""
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_global_confidence.csv")


def extract_global_confidence_dir(
    pae_dir: str,
    *,
    structs_dir: Optional[str] = None,
    save_path: Optional[str] = None,
    ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Global fold confidence (pTM[/ipTM]) for every ``<ID>_pae.npz`` in `pae_dir`.

    Writes a CSV keyed by ID and returns the DataFrame. By default one row per npz in
    `pae_dir`; pass `ids` to score a fixed ID set (missing npz -> NaN row), which lets
    the orchestrator key it off the structs dir like the rest of the structure branch.
    The ``iptm`` column is emitted only if at least one npz carries a non-NaN ipTM
    (single-chain TPS have none). The CSV path defaults to
    ``<structs_dir>_global_confidence.csv`` (or ``<pae_dir>_global_confidence.csv`` if
    no structs dir is given)."""
    if ids is None:
        ids = _ids_from_pae_dir(pae_dir)
    if not ids:
        raise ValueError(
            f"No <ID>_pae.npz files found in {pae_dir} "
            "(written at fold time by esmfold.py / scripts/alphafold/extract_pae.py)."
        )
    print(f"{len(ids)} design(s); reading global confidence from {pae_dir}")

    rows: List[Dict[str, object]] = []
    n_missing = 0
    n_no_ptm = 0
    for design_id in ids:
        path = _pae_path(pae_dir, design_id)
        if not os.path.isfile(path):
            n_missing += 1
        try:
            vals = load_global_confidence(path)
        except Exception as exc:  # malformed npz -> NaN row, keep going
            print(f"  [warn] {design_id}: {exc}")
            vals = {"ptm": float("nan"), IPTM_COLUMN: float("nan")}
        if os.path.isfile(path) and np.isnan(vals["ptm"]):
            n_no_ptm += 1
        rows.append({"ID": design_id, **vals})

    has_iptm = any(not np.isnan(r[IPTM_COLUMN]) for r in rows)
    columns = BASE_COLUMNS + ([IPTM_COLUMN] if has_iptm else [])
    df = pd.DataFrame(rows).reindex(columns=columns)
    df = df.sort_values("ID").reset_index(drop=True)

    print(
        f"Read pTM for {len(df) - n_missing - n_no_ptm} design(s); "
        f"{n_no_ptm} npz without ptm; {n_missing} missing npz (NaN)."
        + (" ipTM present." if has_iptm else " (single-chain: no ipTM column).")
    )

    if save_path is None:
        save_path = (
            _default_save_path(structs_dir) if structs_dir else _default_save_path(pae_dir)
        )
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
