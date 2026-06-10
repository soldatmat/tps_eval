"""TPS structural-domain composition per design, via EnzymeExplorer.

For each input structure this records HOW MANY and WHICH TYPES of TPS
structural domains EnzymeExplorer (EE) detects, as a CSV keyed by ``ID``
(ID == structure filename stem) usable as a filtration criterion alongside
the other tps_eval metrics.

How EE detects/records domains (the key mapping):
  * EE's ``detect_domains`` entry point
    (``enzymeexplorer.src.structure_processing.domain_detections``) aligns each
    structure against seven curated TPS domain TEMPLATES (alpha, beta, gamma,
    ids, delta, epsilon, zeta) with PyMOL + foldseek and returns a mapping
    ``{structure_stem: [MappedRegion, ...]}``, where each ``MappedRegion`` has a
    ``.domain`` (the type) and a ``.module_id`` of the form
    ``"<stem>_<type>_<index>"``. It ALSO writes a portable JSON sidecar
    (``<detections>.json``) via ``save_seq_to_regions_json`` with the same shape
    — no EE import needed to parse it.
  * CRITICAL: the returned mapping is a ``defaultdict(list)`` populated only when
    a region is detected, so structures with ZERO detected domains are simply
    ABSENT from it (and from the JSON). We therefore enumerate the FULL input ID
    set independently (from the structures dir) and left-join, so every design
    gets exactly one row — zero-domain designs included (``n_domains=0``, all
    per-type counts 0, empty ``domain_architecture``).

Detection is PyMOL + foldseek only (no PLM/ESM): CPU-only, no GPU needed.
"""

from __future__ import annotations

import json
import os
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

# The seven canonical TPS structural-domain template types EE detects, in a
# fixed order (matches DEFAULT_DOMAIN_TEMPLATES in EE's domain_detections.py).
# "ids" == the isoprenyl-diphosphate-synthase-like domain.
DOMAIN_TYPES: List[str] = ["alpha", "beta", "gamma", "ids", "delta", "epsilon", "zeta"]

_COUNT_COLS = [f"n_{t}" for t in DOMAIN_TYPES]
COLUMNS = ["ID", "n_domains"] + _COUNT_COLS + ["domain_architecture"]


def _structure_ids(structs_dir: str) -> List[str]:
    """Every input design's ID (== filename stem) in `structs_dir`.

    EE's domain detection consumes ``.pdb`` files (it globs ``*.pdb``); we mirror
    that so the ID universe matches exactly what detection saw. A ``.cif``-only
    structure would not be folded by EE here, so it is not an input to this
    metric. If both .pdb and .cif exist for a stem, the .pdb is the one detected.
    """
    stems = sorted({p.stem for p in Path(structs_dir).glob("*.pdb")})
    return stems


def regions_to_rows(
    seq_to_regions: Dict[str, List[dict]],
    all_ids: List[str],
) -> pd.DataFrame:
    """Build the per-design composition DataFrame.

    `seq_to_regions` maps ``ID -> [region-dict, ...]`` (each region dict has at
    least a ``"domain"`` key and, optionally, an ordering hint). `all_ids` is the
    FULL list of input designs — IDs absent from `seq_to_regions` become
    zero-domain rows. Returns a DataFrame with exactly one row per ID in
    `all_ids`, columns = `COLUMNS`.
    """
    rows: List[Dict[str, object]] = []
    for design_id in all_ids:
        regions = seq_to_regions.get(design_id, []) or []
        # Order regions for the architecture string. EE issues module_ids as
        # "<stem>_<type>_<index>"; sort by (type, index) for a stable, readable
        # ordering. Fall back to insertion order if module_id is missing.
        def _key(r: dict):
            mid = r.get("module_id", "")
            parts = mid.rsplit("_", 1)
            idx = int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else 0
            return (r.get("domain", ""), idx)

        ordered = sorted(regions, key=_key)
        types = [r.get("domain", "") for r in ordered]
        counts = {f"n_{t}": 0 for t in DOMAIN_TYPES}
        for t in types:
            col = f"n_{t}"
            if col in counts:
                counts[col] += 1
            # An unexpected domain type (not in DOMAIN_TYPES) is still counted in
            # n_domains and the architecture string, just without its own column.
        row: Dict[str, object] = {"ID": design_id, "n_domains": len(types)}
        row.update(counts)
        row["domain_architecture"] = "-".join(types)  # "" for zero-domain designs
        rows.append(row)

    df = pd.DataFrame(rows, columns=COLUMNS)
    # Integer dtypes for the count columns (no NaN possible — every ID is filled).
    for col in ["n_domains"] + _COUNT_COLS:
        df[col] = df[col].astype(int)
    return df.sort_values("ID").reset_index(drop=True)


def load_detections_json(json_path: str) -> Dict[str, List[dict]]:
    """Parse EE's portable domain-detection JSON sidecar into ``ID -> [region]``."""
    with open(json_path) as f:
        payload = json.load(f)
    # Shape: {seq_id: [{"module_id":..., "domain":..., "tmscore":..., ...}, ...]}
    return {k: list(v) for k, v in payload.items()}


def detect_domains_json(
    structs_dir: str,
    *,
    out_json_path: str,
    n_jobs: int = 10,
    n_iters: int = 3,
) -> Dict[str, List[dict]]:
    """Run EE domain detection on `structs_dir` and return ``ID -> [region]``.

    Imports EnzymeExplorer (must be installed in the active env) and calls its
    public ``detect_domains`` with an explicit args namespace — this bypasses
    both EE's configargparse default-config requirement (a ``configs/...yaml``
    that need not exist) and the ``run_domain_detection`` wrapper's
    ``str(None)`` bug for ``needed_proteins_csv_path``. Writes EE's pickle +
    JSON sidecar next to `out_json_path`; we parse the JSON.
    """
    import argparse as _argparse

    import yaml  # provided by the EE env
    from enzymeexplorer.src.structure_processing.domain_detections import (
        DEFAULT_DOMAIN_TEMPLATES,
        detect_domains,
    )

    out_json = Path(out_json_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    # EE writes the JSON sidecar as <detections_output_path>.json — so point the
    # pickle at <stem>.pkl so the sidecar lands exactly at out_json_path.
    detections_pkl = out_json.with_suffix(".pkl")
    domains_dir = out_json.parent / "_ee_domains_scratch"
    ssr_pkl = out_json.parent / "_ee_secondary_structure_residues.pkl"

    templates = [yaml.safe_dump(t) for t in DEFAULT_DOMAIN_TEMPLATES]
    args = _argparse.Namespace(
        input_directory_with_structures=str(structs_dir),
        needed_proteins_csv_path=None,  # -> EE globs *.pdb in the input dir
        csv_id_column=None,
        detections_output_path=str(detections_pkl),
        detected_regions_root_path=None,
        domains_output_path=str(domains_dir),
        n_jobs=n_jobs,
        n_iters=n_iters,
        is_bfactor_confidence=True,  # generated structs carry pLDDT in B-factor
        do_not_store_intermediate_files=True,
        store_domains=True,
        detect_multiple_domains_in_each_iteration=True,
        secondary_structure_residues_path=str(ssr_pkl),
        recompute_existing_secondary_structure_residues=True,
        prefilter_pdbs_by_foldseek=False,
        prefilter_e_value=10.0,
        postfilter_domains_by_foldseek=False,
        postfilter_e_value=5.0,
        domain_templates=templates,
    )
    detect_domains(args)
    return load_detections_json(str(out_json))


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_domain_composition.csv")


def extract_domain_composition_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    detections_json: Optional[str] = None,
    n_jobs: int = 10,
    n_iters: int = 3,
) -> pd.DataFrame:
    """Per-design TPS-domain composition for every structure in `structs_dir`.

    If `detections_json` is given AND exists, parse that EE domain-detection JSON
    sidecar instead of re-running detection (cheap re-use of a prior run).
    Otherwise run EE's ``detect_domains`` to produce one. Either way, EVERY input
    design (every ``*.pdb`` stem) gets exactly one row — including designs with
    zero detected domains. Writes the CSV keyed by ID and returns the DataFrame.
    """
    all_ids = _structure_ids(structs_dir)
    if not all_ids:
        raise ValueError(
            f"No .pdb structures found in {structs_dir} "
            "(EE domain detection consumes .pdb files; ID = filename stem)."
        )
    print(f"{len(all_ids)} input design(s) in {structs_dir}")

    if detections_json and os.path.isfile(detections_json):
        print(f"Reusing existing domain detections from {detections_json}")
        seq_to_regions = load_detections_json(detections_json)
    else:
        if save_path is None:
            save_path = _default_save_path(structs_dir)
        json_path = detections_json or (
            os.path.splitext(save_path)[0] + "_detections.json"
        )
        print(f"Running EnzymeExplorer domain detection -> {json_path}")
        seq_to_regions = detect_domains_json(
            structs_dir, out_json_path=json_path, n_jobs=n_jobs, n_iters=n_iters
        )

    n_with = sum(1 for i in all_ids if seq_to_regions.get(i))
    print(
        f"EE detected domains in {n_with}/{len(all_ids)} design(s); "
        f"{len(all_ids) - n_with} design(s) have zero domains (kept as n_domains=0 rows)."
    )

    df = regions_to_rows(seq_to_regions, all_ids)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
