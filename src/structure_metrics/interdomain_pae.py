"""Inter-domain PAE: confidence in the RELATIVE orientation of a design's domains.

The Predicted Aligned Error (PAE) matrix from AlphaFold/ESMFold reports, for every
ordered residue pair (i, j), the expected position error (in Angstrom) of residue j
when the predicted and true structures are superposed on residue i. The DIAGONAL
blocks (residues within one domain) are usually confident even for nonsense
multi-domain designs — each domain folds locally. The OFF-DIAGONAL blocks (residue
in domain A vs residue in domain B) report whether the model is confident about how
the two domains are PACKED relative to each other. A high inter-domain PAE with low
pLDDT-within-domain is the classic "two confident domains floating in uncertain
relative orientation" failure mode that per-residue pLDDT cannot see.

This tool, per design:
  1. Gets per-domain RESIDUE RANGES from the EnzymeExplorer (EE) domain detector —
     the SAME ``detect_domains`` call ``domain_composition.py`` uses, but it KEEPS the
     residue spans (the integer keys of each region's ``residues_mapping``, which are
     PDB author residue numbers) rather than collapsing to per-type counts.
  2. Loads ``<ID>_pae.npz`` (written at fold time by esmfold.py / the AF3 extractor),
     maps each domain's residue numbers onto PAE matrix indices via the npz's
     ``residue_ids`` axis, and reduces the inter-domain blocks: for every pair of
     domains (A, B) it averages PAE over residues(A) x residues(B) in BOTH directions
     (PAE is asymmetric) and averages the two means.

Per design (keyed by ``ID``) it emits:
  * ``mean_interdomain_pae`` — mean over all inter-domain domain-pairs.
  * ``max_interdomain_pae``  — worst (largest) inter-domain pair.
  * ``n_domains``            — number of EE-detected domains used.
  * optional per-pair columns ``pae_<A>_<B>`` (use ``--per_pair``).

A SINGLE-domain (or zero-domain) design has no inter-domain block: it is emitted as
N/A (NaN, ``n_domains`` recorded). NaN is also emitted when the PAE npz is missing or
domains could not be detected. Raw numbers only — "natural TPS" bands are computed
separately by the reference-stats pipeline.

The PAE npz schema (one schema, both folders — see esmfold.py / the AF3 extractor):
  * ``pae``         : float32 (L, L), Angstrom; PAE[i, j] aligned-on-i error of j.
  * ``residue_ids`` : int32 (L,), the PDB author residue number of each PAE row/col.
  * ``n_residues``  : int scalar (== L).
  * ``source``      : str, 'esmfold' | 'alphafold3'.
"""

from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Reuse the canonical structure-loading / ID-stem conventions and EE domain
# detection (residue spans) so this tool stays source-agnostic and mirrors the
# rest of the structure branch.
from enzyme_explorer.domain_composition import (  # noqa: E402
    _structure_ids,
    detect_domains_json,
    load_detections_json,
)

# Base (always present) output columns; per-pair columns are appended on demand.
BASE_COLUMNS = ["ID", "mean_interdomain_pae", "max_interdomain_pae", "n_domains"]


# --------------------------------------------------------------------------- #
# PAE npz loading                                                             #
# --------------------------------------------------------------------------- #
def _pae_path(pae_dir: str, design_id: str) -> str:
    return os.path.join(pae_dir, f"{design_id}_pae.npz")


def load_pae(pae_path: str):
    """Load a ``<ID>_pae.npz`` -> (pae (L,L) float array, residue_ids (L,) int array).

    Returns ``(None, None)`` if the file is missing. Raises on a malformed file
    (shape mismatch) so a corrupt artefact is loud rather than silently NaN.
    """
    if not os.path.isfile(pae_path):
        return None, None
    with np.load(pae_path, allow_pickle=False) as npz:
        pae = np.asarray(npz["pae"], dtype=float)
        if "residue_ids" in npz.files:
            residue_ids = np.asarray(npz["residue_ids"]).astype(int)
        else:
            # Older/minimal npz without an explicit axis: assume PAE rows are the
            # contiguous residues 1..L in PDB author numbering (the ESMFold case).
            residue_ids = np.arange(1, pae.shape[0] + 1, dtype=int)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE in {pae_path} is not a square 2D matrix (shape {pae.shape})")
    if residue_ids.shape[0] != pae.shape[0]:
        raise ValueError(
            f"residue_ids length {residue_ids.shape[0]} != PAE dimension {pae.shape[0]} "
            f"in {pae_path}"
        )
    return pae, residue_ids


# --------------------------------------------------------------------------- #
# Domain residue ranges -> PAE indices                                        #
# --------------------------------------------------------------------------- #
def domain_residue_numbers(regions: List[dict]) -> Dict[str, List[int]]:
    """Map each EE region's ``module_id`` -> sorted list of PDB residue numbers.

    The residue span of a region is the set of integer KEYS of its
    ``residues_mapping`` (PyMOL ``resi`` values == PDB author residue numbers).
    We key by ``module_id`` (unique "<stem>_<type>_<index>") so two domains of the
    same type ("alpha", "alpha") stay distinct domains.
    """
    out: Dict[str, List[int]] = {}
    for region in regions:
        mapping = region.get("residues_mapping", {}) or {}
        resis = sorted({int(k) for k in mapping.keys()})
        if not resis:
            continue
        module_id = region.get("module_id") or f"{region.get('domain', 'dom')}_{len(out)}"
        out[module_id] = resis
    return out


def _resis_to_indices(resis: List[int], id_to_index: Dict[int, int]) -> List[int]:
    """Map PDB residue numbers onto PAE matrix indices via the npz residue_ids axis,
    dropping any residue absent from the PAE axis (defensive against minor numbering
    mismatches between EE detection and the folded structure)."""
    return [id_to_index[r] for r in resis if r in id_to_index]


def interdomain_pae_blocks(
    pae: np.ndarray,
    residue_ids: np.ndarray,
    domain_to_resis: Dict[str, List[int]],
) -> Dict[str, float]:
    """Mean PAE for every ordered domain pair, symmetrised.

    For domains A, B the inter-domain score is
        0.5 * (mean(PAE[idx(A), idx(B)]) + mean(PAE[idx(B), idx(A)]))
    i.e. the off-diagonal A-B and B-A blocks averaged (PAE is asymmetric).
    Returns ``{"<A>_<B>": score}`` for every unordered pair with non-empty index
    sets on both sides. Empty when fewer than 2 domains map onto the matrix.
    """
    id_to_index = {int(r): i for i, r in enumerate(residue_ids.tolist())}
    dom_indices = {
        d: _resis_to_indices(resis, id_to_index) for d, resis in domain_to_resis.items()
    }
    # Keep only domains that actually land on the PAE axis.
    dom_indices = {d: idx for d, idx in dom_indices.items() if idx}

    pairs: Dict[str, float] = {}
    for a, b in combinations(sorted(dom_indices), 2):
        ia = np.asarray(dom_indices[a], dtype=int)
        ib = np.asarray(dom_indices[b], dtype=int)
        block_ab = pae[np.ix_(ia, ib)]  # aligned-on-A error of B residues
        block_ba = pae[np.ix_(ib, ia)]  # aligned-on-B error of A residues
        pairs[f"{a}_{b}"] = 0.5 * (float(block_ab.mean()) + float(block_ba.mean()))
    return pairs


# --------------------------------------------------------------------------- #
# Per-design row + directory driver                                           #
# --------------------------------------------------------------------------- #
def _na_row(design_id: str, n_domains: Optional[int]) -> Dict[str, object]:
    return {
        "ID": design_id,
        "mean_interdomain_pae": np.nan,
        "max_interdomain_pae": np.nan,
        "n_domains": (np.nan if n_domains is None else int(n_domains)),
    }


def design_row(
    design_id: str,
    pae_dir: str,
    regions: List[dict],
    *,
    per_pair: bool = False,
) -> Dict[str, object]:
    """One output row for a single design.

    N/A (NaN) when: the PAE npz is missing; EE detected fewer than 2 domains; or no
    pair of domains maps onto the PAE matrix. ``n_domains`` is recorded whenever
    domains were detected (even single-domain N/A rows), and is NaN only when the
    PAE artefact is missing (we then don't know the residue count to trust).
    """
    domain_to_resis = domain_residue_numbers(regions)
    n_domains = len(domain_to_resis)

    pae, residue_ids = load_pae(_pae_path(pae_dir, design_id))
    if pae is None:
        # No PAE artefact -> nothing computable. Still record n_domains (EE ran) so a
        # single-domain vs multi-domain-but-missing-PAE N/A can be told apart.
        return _na_row(design_id, n_domains)

    # Single- or zero-domain designs have no inter-domain block by definition.
    if n_domains < 2:
        return _na_row(design_id, n_domains)

    pairs = interdomain_pae_blocks(pae, residue_ids, domain_to_resis)
    if not pairs:
        return _na_row(design_id, n_domains)

    values = np.array(list(pairs.values()), dtype=float)
    row: Dict[str, object] = {
        "ID": design_id,
        "mean_interdomain_pae": float(values.mean()),
        "max_interdomain_pae": float(values.max()),
        "n_domains": int(n_domains),
    }
    if per_pair:
        for pair_name, score in pairs.items():
            row[f"pae_{pair_name}"] = score
    return row


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_interdomain_pae.csv")


def extract_interdomain_pae_dir(
    structs_dir: str,
    pae_dir: str,
    *,
    save_path: Optional[str] = None,
    detections_json: Optional[str] = None,
    per_pair: bool = False,
    n_jobs: int = 10,
    n_iters: int = 3,
) -> pd.DataFrame:
    """Inter-domain PAE for every design in `structs_dir`, paired with its PAE npz in
    `pae_dir`. Writes a CSV keyed by ID and returns the DataFrame. EVERY input design
    (every ``*.pdb`` stem) gets exactly one row; single-domain / missing-PAE /
    no-domain designs are N/A (NaN)."""
    all_ids = _structure_ids(structs_dir)
    if not all_ids:
        raise ValueError(
            f"No .pdb structures found in {structs_dir} "
            "(EE domain detection consumes .pdb files; ID = filename stem)."
        )
    print(f"{len(all_ids)} input design(s) in {structs_dir}; PAE dir {pae_dir}")

    # Domain residue ranges via EE detect_domains — reuse domain_composition's call.
    if detections_json and os.path.isfile(detections_json):
        print(f"Reusing existing domain detections from {detections_json}")
        seq_to_regions = load_detections_json(detections_json)
    else:
        if save_path is None:
            save_path = _default_save_path(structs_dir)
        json_path = detections_json or (os.path.splitext(save_path)[0] + "_detections.json")
        print(f"Running EnzymeExplorer domain detection -> {json_path}")
        seq_to_regions = detect_domains_json(
            structs_dir, out_json_path=json_path, n_jobs=n_jobs, n_iters=n_iters
        )

    rows: List[Dict[str, object]] = []
    n_missing_pae = 0
    n_single = 0
    n_scored = 0
    for design_id in all_ids:
        regions = seq_to_regions.get(design_id, []) or []
        try:
            row = design_row(design_id, pae_dir, regions, per_pair=per_pair)
        except Exception as exc:  # malformed npz -> NaN row, keep going
            print(f"  [warn] {design_id}: {exc}")
            row = _na_row(design_id, len(domain_residue_numbers(regions)))
        rows.append(row)
        if not os.path.isfile(_pae_path(pae_dir, design_id)):
            n_missing_pae += 1
        elif pd.isna(row["mean_interdomain_pae"]):
            n_single += 1
        else:
            n_scored += 1

    # Stable column order: base columns first, then any per-pair columns sorted.
    extra = sorted({c for r in rows for c in r if c not in BASE_COLUMNS})
    columns = BASE_COLUMNS + extra
    df = pd.DataFrame(rows).reindex(columns=columns)
    df["n_domains"] = df["n_domains"].astype("Int64")  # nullable int (NaN for missing PAE)
    df = df.sort_values("ID").reset_index(drop=True)

    print(
        f"Scored inter-domain PAE for {n_scored} design(s); "
        f"{n_single} single/zero-domain N/A; {n_missing_pae} missing PAE npz."
    )
    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
