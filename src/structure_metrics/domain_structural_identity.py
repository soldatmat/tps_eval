"""Domain-level structural-identity of each generated design to known TPS domains.

This is the *domain-level* analog of the full-structure ``structural_identity``
metric (``src/structure_metrics/run_structural_identity.py``). Instead of aligning
the whole generated structure against whole known-TPS structures, it:

  1. Detects the individual TPS structural domains in each generated structure
     (via EnzymeExplorer's ``detect_domains``), writing one ``.pdb`` per detected
     domain.
  2. Foldseek-aligns those detected domain structures against EnzymeExplorer's
     curated reference set of *known-TPS* domain structures (the martsDB detected
     domains), via ``src/foldseek/domain_alignment.py``.
  3. Reduces to one row per design (keyed by ``ID``): the best domain-level
     TM-score (+ lddt) to the nearest known-TPS domain, which reference domain it
     matched, that reference's domain TYPE, and per-domain-type bests.

A design with NO detected domain gets a single NaN row (kept, so the ID universe
matches the rest of the structure branch).

How EE detect_domains output feeds domain_alignment (the key plumbing)
----------------------------------------------------------------------
EE's ``detect_domains`` (``store_domains=True``) writes each detected domain as a
``.pdb`` named ``<ID>_<type>_<index>.pdb`` (``module_id``) — e.g.
``mydesign_alpha_0.pdb`` — under the ``domains_output_path`` we give it, in BOTH a
flat layout (``<domains_output_path>/<module_id>.pdb``) and a per-type layout
(``<domains_output_path>/<type>/<module_id>.pdb``). So the domain TYPE is encoded
in the FILENAME, both for the detected domains and for the reference domains
(``marts_E00000_alpha_0.pdb`` ...). We therefore:

  * feed the FLAT detected-domains dir as ``domain_alignment``'s
    ``--detected_domain_structures_root``;
  * feed the FLAT EE reference-domains directory (which holds all ~2.4k
    ``marts_*_<type>_*.pdb`` files) as ``--known_domain_structures_root``;
  * parse the foldseek query stem (``<ID>_<type>_<index>``) back to the design ID,
    and the target stem (``marts_<acc>_<type>_<index>``) to the reference's TYPE.

This mirrors ``enzyme_explorer/domain_composition.py``'s use of the SAME installed
``detect_domains`` entry point (template base PDBs are auto-resolved by EE from its
own ``data/domain_templates``; the per-domain save runs inside EE's spawn pool with
cwd == the input directory, so no template staging or cwd juggling is needed here).
The portable detections JSON sidecar (``<detections>.json``) gives the per-design
detected-domain COUNT so designs with zero domains become NaN rows.
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from foldseek.domain_alignment import main as _domain_alignment  # noqa: E402

# The TPS structural-domain template types, fixed order (matches
# enzyme_explorer/domain_composition.DOMAIN_TYPES). The TYPE is parsed from the
# detected/reference domain filenames "<stem>_<type>_<index>".
DOMAIN_TYPES: List[str] = ["alpha", "beta", "gamma", "ids", "delta", "epsilon", "zeta"]

# Per-design output columns. The "to_known" naming mirrors structural_identity.
_BASE_COLUMNS = [
    "ID",
    "domain_structural_tmscore_to_known",
    "domain_structural_tmscore_to_known_hit",
    "domain_structural_tmscore_to_known_type",
    "domain_structural_lddt_to_known",
    "n_detected_domains",
]
_PER_TYPE_COLUMNS = [f"domain_structural_tmscore_to_known_{t}" for t in DOMAIN_TYPES]
COLUMNS = _BASE_COLUMNS + _PER_TYPE_COLUMNS


def _stem(name: str) -> str:
    """Filename stem: drop any path and a trailing structure extension."""
    base = os.path.basename(str(name))
    for ext in (".pdb.gz", ".cif.gz", ".pdb", ".cif", ".ent"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def _parse_module_id(stem: str) -> tuple[str, str]:
    """Split a domain ``module_id`` stem ``<id>_<type>_<index>`` into (id, type).

    The id may itself contain underscores (designs are often ``run_0001`` etc.),
    so we anchor on the *known* domain TYPE token that EE always inserts as the
    second-to-last underscore-field: ``<id>_<type>_<index>``. Falls back to
    ("<stem>", "") if the pattern is not recognised.
    """
    parts = stem.rsplit("_", 2)
    if len(parts) == 3 and parts[2].isdigit() and parts[1] in DOMAIN_TYPES:
        return parts[0], parts[1]
    # Type token not recognised — treat the whole stem as the id (no type).
    return stem, ""


def _structure_ids(structs_dir: str) -> List[str]:
    """Every input design's ID (== .pdb filename stem) in `structs_dir`. EE domain
    detection consumes ``.pdb`` files; mirror enzyme_explorer/domain_composition so
    the ID universe matches exactly what detection saw."""
    return sorted({_stem(p) for p in glob.glob(os.path.join(structs_dir, "*.pdb"))})


def detect_domains_to_dir(
    structs_dir: str,
    *,
    out_domains_dir: str,
    n_jobs: int = 10,
    n_iters: int = 3,
) -> Dict[str, int]:
    """Run EE ``detect_domains`` on `structs_dir`, writing per-domain ``.pdb``
    structures into `out_domains_dir` (flat ``<module_id>.pdb`` + per-type subdirs).

    Returns ``ID -> n_detected_domains`` for every input design (zero-domain
    designs included, value 0), read from EE's portable detections JSON sidecar.
    Imports EnzymeExplorer (must be installed in the active env, alongside
    foldseek). Uses the installed ``detect_domains`` entry point exactly as
    ``enzyme_explorer/domain_composition.py`` does — templates are auto-resolved by
    EE, and the per-domain save runs inside EE's spawn pool (cwd == input dir).
    """
    import argparse as _argparse

    import yaml  # provided by the EE env
    from enzymeexplorer.src.structure_processing.domain_detections import (
        DEFAULT_DOMAIN_TEMPLATES,
        detect_domains,
    )

    all_ids = _structure_ids(structs_dir)
    if not all_ids:
        raise ValueError(
            f"No .pdb structures found in {structs_dir} "
            "(EE domain detection consumes .pdb files; ID = filename stem)."
        )

    out_domains = Path(out_domains_dir).absolute()
    out_domains.mkdir(parents=True, exist_ok=True)

    # Scratch for EE's own sidecars (pickle/json/ss-cache). The detected per-domain
    # pdbs go to out_domains (absolute, so EE's `Path(cwd)/domains_output_path` join
    # resolves to it unchanged).
    work = Path(tempfile.mkdtemp(prefix="domain_struct_id_"))
    detections_pkl = work / "detections.pkl"
    ssr_pkl = work / "_ssr.pkl"

    templates = [yaml.safe_dump(t) for t in DEFAULT_DOMAIN_TEMPLATES]
    args = _argparse.Namespace(
        input_directory_with_structures=str(structs_dir),
        needed_proteins_csv_path=None,  # -> EE globs *.pdb in the input dir
        csv_id_column=None,
        detections_output_path=str(detections_pkl),
        detected_regions_root_path=None,
        domains_output_path=str(out_domains),
        n_jobs=n_jobs,
        n_iters=n_iters,
        is_bfactor_confidence=True,  # generated structs carry pLDDT in B-factor
        do_not_store_intermediate_files=True,
        store_domains=True,  # WRITE the per-domain .pdb structures
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

    # Per-design detected-domain counts from the portable JSON sidecar (so designs
    # with zero detected domains — ABSENT from the mapping — still get a row).
    counts: Dict[str, int] = {sid: 0 for sid in all_ids}
    json_path = detections_pkl.with_suffix(".json")
    if json_path.is_file():
        import json
        with open(json_path) as f:
            payload = json.load(f)
        for seq_id, regions in payload.items():
            counts[seq_id] = counts.get(seq_id, 0) + len(regions or [])

    n_flat = len(glob.glob(os.path.join(str(out_domains), "*.pdb")))
    print(
        f"EE detected {n_flat} domain structure(s) across "
        f"{sum(1 for v in counts.values() if v)}/{len(all_ids)} design(s) "
        f"-> {out_domains}"
    )
    shutil.rmtree(work, ignore_errors=True)
    return counts


def _reduce_alignments(
    raw_hits_csv: str,
    domain_counts: Dict[str, int],
    *,
    exclude_self: bool = False,
) -> pd.DataFrame:
    """Reduce per-detected-domain foldseek hits to one row per design ID.

    `raw_hits_csv` is the full per-hit table from ``domain_alignment`` (one row per
    query-domain × reference-domain hit). `domain_counts` maps every input ID to
    its number of detected domains (so zero-domain designs become NaN rows).

    `exclude_self`: searching a domain set against itself (leave-one-out). Drop hits
    to a reference domain originating from the query's OWN source structure (the
    reference domain stem ``<id>_<type>_<index>`` parses to the same design id),
    BEFORE the best-hit reduction, so each design's best hit is its nearest OTHER
    known-TPS domain instead of the trivial self-match TM~1.0.
    """
    hits = pd.read_csv(raw_hits_csv)

    # Map each hit to (design_id, ref_type) from the query/target filenames.
    hits = hits.copy()
    hits["__query_stem"] = hits["query"].map(_stem)
    hits["__target_stem"] = hits["target"].map(_stem)
    hits["__design_id"] = hits["__query_stem"].map(lambda s: _parse_module_id(s)[0])
    hits["__ref_type"] = hits["__target_stem"].map(lambda s: _parse_module_id(s)[1])
    if exclude_self:
        # "self" == the reference domain comes from this design's own source
        # structure (same source id as the query domain's design id).
        hits["__ref_source_id"] = hits["__target_stem"].map(
            lambda s: _parse_module_id(s)[0]
        )
        hits = hits[hits["__ref_source_id"] != hits["__design_id"]]
    # alntmscore is the symmetric domain-level TM-score (LARGER = closer).
    hits["alntmscore"] = pd.to_numeric(hits["alntmscore"], errors="coerce")
    hits["lddt"] = pd.to_numeric(hits["lddt"], errors="coerce")

    rows: List[Dict[str, object]] = []
    by_design = {d: g for d, g in hits.groupby("__design_id", sort=False)}
    for design_id in sorted(domain_counts):
        row: Dict[str, object] = {col: np.nan for col in COLUMNS}
        row["ID"] = design_id
        row["n_detected_domains"] = int(domain_counts[design_id])
        grp = by_design.get(design_id)
        if grp is not None and len(grp) and grp["alntmscore"].notna().any():
            best_idx = grp["alntmscore"].idxmax()
            best = grp.loc[best_idx]
            row["domain_structural_tmscore_to_known"] = float(best["alntmscore"])
            row["domain_structural_tmscore_to_known_hit"] = str(best["__target_stem"])
            row["domain_structural_tmscore_to_known_type"] = str(best["__ref_type"])
            row["domain_structural_lddt_to_known"] = (
                float(best["lddt"]) if pd.notna(best["lddt"]) else np.nan
            )
            # Per-domain-type best (by the reference domain's type).
            for ref_type, tg in grp.groupby("__ref_type", sort=False):
                col = f"domain_structural_tmscore_to_known_{ref_type}"
                if col in row and tg["alntmscore"].notna().any():
                    row[col] = float(tg["alntmscore"].max())
        rows.append(row)

    df = pd.DataFrame(rows, columns=COLUMNS)
    df["n_detected_domains"] = df["n_detected_domains"].astype(int)
    return df.sort_values("ID").reset_index(drop=True)


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(
        os.path.dirname(d), os.path.basename(d) + "_domain_structural_identity.csv"
    )


def extract_domain_structural_identity_dir(
    structs_dir: str,
    known_domain_structures_root: str,
    *,
    save_path: Optional[str] = None,
    n_jobs: int = 10,
    n_iters: int = 3,
    keep_detected_domains: Optional[str] = None,
    exclude_self: bool = False,
) -> pd.DataFrame:
    """Domain-level structural identity for every design in `structs_dir`.

    Detects TPS domains, foldseek-aligns them against the known-TPS reference
    domains under `known_domain_structures_root`, reduces to one row per design ID,
    writes ``<structs_dir>_domain_structural_identity.csv`` and returns the df.

    `exclude_self`: self-search (leave-one-out). Before the best-hit reduction, drop
    hits to a reference domain originating from the query's own source structure, so
    a domain set searched against itself yields the nearest OTHER known-TPS domain
    instead of the trivial self-match TM~1.0. Defaults OFF (gen-vs-reference runs).
    """
    if not os.path.isdir(known_domain_structures_root):
        raise NotADirectoryError(
            f"known_domain_structures_root does not exist: {known_domain_structures_root}"
        )

    detected_root = keep_detected_domains
    tmp_detected: Optional[str] = None
    if detected_root is None:
        tmp_detected = tempfile.mkdtemp(prefix="detected_domains_")
        detected_root = tmp_detected
    else:
        os.makedirs(detected_root, exist_ok=True)

    tmp_align = tempfile.mkdtemp(prefix="domain_alignment_")
    # EE writes each detected domain BOTH flat (<root>/<module_id>.pdb) AND under a
    # per-type subdir (<root>/<type>/<module_id>.pdb). foldseek easy-search recurses
    # into subdirs, so to avoid double-counting we align against a FLAT-ONLY dir
    # holding just the top-level files (one per detected domain).
    flat_align = tempfile.mkdtemp(prefix="detected_flat_")
    try:
        domain_counts = detect_domains_to_dir(
            structs_dir,
            out_domains_dir=detected_root,
            n_jobs=n_jobs,
            n_iters=n_iters,
        )

        flat_files = sorted(glob.glob(os.path.join(detected_root, "*.pdb")))
        for src in flat_files:
            link = os.path.join(flat_align, os.path.basename(src))
            try:
                os.symlink(os.path.abspath(src), link)
            except OSError:
                shutil.copy(src, link)

        if not flat_files:
            # No domain detected in ANY design -> all-NaN rows, skip foldseek.
            print("No domains detected in any design; all rows will be NaN.")
            rows = [
                {**{c: np.nan for c in COLUMNS}, "ID": sid, "n_detected_domains": 0}
                for sid in sorted(domain_counts)
            ]
            df = pd.DataFrame(rows, columns=COLUMNS)
            df["n_detected_domains"] = df["n_detected_domains"].astype(int)
            df = df.sort_values("ID").reset_index(drop=True)
        else:
            _domain_alignment(
                argparse.Namespace(
                    detected_domain_structures_root=flat_align,
                    known_domain_structures_root=known_domain_structures_root,
                    output_root=tmp_align,
                    store_intermediate_results=True,  # keep raw per-hit table
                    random_run_id=False,
                )
            )
            raw_hits = os.path.join(tmp_align, "domain_alignments.csv")
            df = _reduce_alignments(raw_hits, domain_counts, exclude_self=exclude_self)
    finally:
        shutil.rmtree(tmp_align, ignore_errors=True)
        shutil.rmtree(flat_align, ignore_errors=True)
        if tmp_detected is not None:
            shutil.rmtree(tmp_detected, ignore_errors=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    n_with = int((df["n_detected_domains"] > 0).sum())
    print(
        f"Wrote {len(df)} rows to {save_path} "
        f"({n_with} design(s) with >=1 detected domain)"
    )
    return df
