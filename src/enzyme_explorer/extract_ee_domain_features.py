"""Extract EnzymeExplorer (EE) domain-COMPARISON features — the structure/function
feature block of EE's production model
``PlmDomainsRandomForest__tps_esm-1v-subseq_..._domains_subset``.

This is **block (A)** of the two EE feature blocks (the PLM block is produced by the
sibling ``extract_ee_esm1v_embeddings.py``). It reproduces, for each input protein, the
exact ``dom_feat`` vector that EE's ``scripts/easy_predict.py`` feeds to the production
classifier:

  1. Detect TPS structural domains on each AF/ESMFold structure
     (``enzymeexplorer.src.structure_processing.domain_detections``).
  2. Compare every detected domain to EE's curated reference functional-domain modules
     with Foldseek TM-score
     (``enzymeexplorer.src.structure_processing.comparing_to_known_domains_foldseek``),
     restricted to the production ``domains_subset.pkl`` (the ``domains_subset`` the
     production model uses).
  3. For each reference module ``m`` (of type ``t``) the feature value is
     ``1 - max TM-score`` over the protein's detected domains of type ``t`` against ``m``
     (0 where no comparison exists) — identical to ``easy_predict.py``'s
     ``dom_feat = 1 - dom_feat`` after filling from ``comparison_results``.

Output: one row per protein, columns = ``id`` + one column per reference module
(``<known_module_id>``), keyed by ``Enzyme_marts_ID``. The column universe is the UNION
of the reference modules across all 5 production fold-classifiers' selected
``domain_type_2_order_of_domain_modules`` — which equals the 897-column
``feat_idx`` selection inside ``domains_subset.pkl``. Fold-specific column subsetting is
left to the classifier; this matrix is the full, deterministic structure/function block.

CRITICAL — the underscore-in-ID bug. EE keys its comparison results on
``row['query'].split('_')[0]`` (``comparing_to_known_domains_foldseek.py``). Detected-
domain PDB filenames are ``{protein_id}_{domain}_{i}.pdb``, so an input ``protein_id``
containing ``_`` (e.g. ``marts_E00000``) collapses to the prefix (``marts``) and EVERY
protein loses its comparisons. We therefore run detection/comparison on **sanitized,
underscore-free IDs** (``marts_E00000`` -> ``martsE00000``) and map the results back to
the original ``Enzyme_marts_ID`` for the output CSV.

Must run in the EE ``enzyme_explorer`` env (PyMOL + Foldseek). CPU-only — no GPU.
Run from the EE ``scripts/`` dir so ``data/`` resolves to the production
``scripts/data/`` bundle (reference domains, classifiers, ``domains_subset.pkl``).
"""

from __future__ import annotations

import argparse
import os
import pickle
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Canonical domain-type order for a stable column layout. EE's classifier renames the
# first/second detected alpha domain to alpha1/alpha2; reference modules are typed
# 'alpha' in the PDB filename but bucketed alpha1/alpha2 in the classifier layout.
DOMAIN_TYPE_ORDER = ["alpha1", "alpha2", "beta", "gamma", "delta", "epsilon"]

_SAFE_RE = re.compile(r"^[A-Za-z0-9]+$")


def _sanitize_id(orig_id: str) -> str:
    """Strip everything that is not [A-Za-z0-9] so the detected-domain filename
    ``{id}_{domain}_{i}`` survives EE's ``split('_')[0]`` keying intact."""
    return re.sub(r"[^A-Za-z0-9]", "", orig_id)


def build_module_columns(
    fold_classifiers: list,
) -> tuple[list[str], dict[str, str]]:
    """Canonical ordered list of reference-module feature columns = the UNION of the
    selected modules across all fold classifiers (== the 897-col ``domains_subset``
    selection). Returns (ordered_module_ids, module_id -> domain_type)."""
    type_to_modules: dict[str, set[str]] = defaultdict(set)
    module_to_type: dict[str, str] = {}
    for clf in fold_classifiers:
        for dtype, module_list in clf.domain_type_2_order_of_domain_modules.items():
            for module_id, _idx in module_list:
                type_to_modules[dtype].add(module_id)
                module_to_type[module_id] = dtype
    ordered: list[str] = []
    for dtype in DOMAIN_TYPE_ORDER:
        ordered.extend(sorted(type_to_modules.get(dtype, set())))
    # any unexpected type (defensive) appended last, sorted by type then module
    for dtype in sorted(set(type_to_modules) - set(DOMAIN_TYPE_ORDER)):
        ordered.extend(sorted(type_to_modules[dtype]))
    return ordered, module_to_type


def stage_structures(
    structs_dir: str, ids: list[str], scratch_dir: Path, *,
    domain_templates_zip: str = "data/domain_templates.zip",
) -> tuple[Path, dict[str, str], pd.DataFrame]:
    """Copy each protein's ``<id>.pdb`` into a scratch dir under a SANITIZED filename.
    Also unzips the domain-template PDBs into the dir — EE's ``domain_detections``
    aligns against them and FileNotFound-errors if they are absent (easy_predict.py
    likewise unzips them into the input dir). The template stems (1ps1, 5eat, 3p5r,
    P48449, Q7Z859) are underscore-free and not in our ``Enzyme_marts_ID`` set, so they
    never enter the output. Returns (sanitized_structs_dir, sanitized_id -> original_id,
    needed_proteins_df)."""
    sdir = scratch_dir / "structs_sanitized"
    sdir.mkdir(parents=True, exist_ok=True)
    if os.path.isfile(domain_templates_zip):
        os.system(f"unzip -o -q {domain_templates_zip} -d {sdir}")
    san_to_orig: dict[str, str] = {}
    rows = []
    src_root = Path(structs_dir)
    for orig in ids:
        src = src_root / f"{orig}.pdb"
        if not src.is_file():
            continue
        san = _sanitize_id(orig)
        # disambiguate accidental collisions
        base = san
        k = 0
        while san in san_to_orig and san_to_orig[san] != orig:
            k += 1
            san = f"{base}x{k}"
        san_to_orig[san] = orig
        shutil.copyfile(src, sdir / f"{san}.pdb")
        rows.append({"ID": san})
    needed_df = pd.DataFrame(rows)
    return sdir, san_to_orig, needed_df


def run_detection_and_comparison(
    san_structs_dir: Path,
    needed_csv: Path,
    scratch_dir: Path,
    *,
    n_jobs: int,
    reuse_cached: bool = False,
) -> tuple[dict, dict]:
    """Run EE domain detection + comparison-to-known-domains on the sanitized structures.
    Returns (detected_domains, comparison_results) — both keyed by SANITIZED id.
    Mirrors the exact CLI calls easy_predict.py makes.

    Pickle paths are DETERMINISTIC (fixed names under ``scratch/_temp``) so a later run
    with ``reuse_cached=True`` can skip the expensive detection+comparison and just
    rebuild the feature matrix (e.g. after a logic fix). The detected-domain region PDBs
    live in a fixed dir likewise.
    """
    temp = scratch_dir / "_temp"
    temp.mkdir(parents=True, exist_ok=True)
    detections_path = temp / "detections.pkl"
    comparison_path = temp / "comparisons.pkl"
    detected_regions_root = temp / "detected_domains"

    if reuse_cached and detections_path.is_file() and comparison_path.is_file():
        print(f"[reuse] loading cached detection+comparison from {temp}", flush=True)
        with open(detections_path, "rb") as f:
            detected_domains = pickle.load(f)
        with open(comparison_path, "rb") as f:
            comparison_results = pickle.load(f)
        return detected_domains, comparison_results

    detected_regions_root.mkdir(parents=True, exist_ok=True)
    py = sys.executable  # use THIS env's python for the inner module calls (not bare `python`)
    det_cmd = (
        f"{py} -m enzymeexplorer.src.structure_processing.domain_detections "
        f'--needed-proteins-csv-path "{needed_csv}" '
        "--csv-id-column ID "
        f"--n-jobs {n_jobs} "
        f"--input-directory-with-structures {san_structs_dir} "
        "--is-bfactor-confidence "
        f'--detections-output-path "{detections_path}" '
        f'--detected-regions-root-path "{detected_regions_root}" '
        f'--domains-output-path "{detected_regions_root}" '
        "--store-domains "
        "--recompute-existing-secondary-structure-residues "
        "--do-not-store-intermediate-files"
    )
    print(f"[detect] {det_cmd}", flush=True)
    rc = os.system(det_cmd)
    if rc != 0:
        raise RuntimeError(f"domain_detections exited with code {rc}")

    with open(detections_path, "rb") as f:
        detected_domains = pickle.load(f)
    print(f"[detect] detected domains in {len(detected_domains)} protein(s)", flush=True)

    comparison_results: dict = {}
    if detected_domains:
        cmp_cmd = (
            f"{py} -m enzymeexplorer.src.structure_processing.comparing_to_known_domains_foldseek "
            "--known-domain-structures-root data/tps_detected_domains/all "
            f'--detected-domain-structures-root "{detected_regions_root}" '
            "--path-to-known-domains-subset data/domains_subset.pkl "
            f'--output-path "{comparison_path}"'
        )
        print(f"[compare] {cmp_cmd}", flush=True)
        rc = os.system(cmp_cmd)
        if rc != 0:
            raise RuntimeError(f"comparing_to_known_domains exited with code {rc}")
        with open(comparison_path, "rb") as f:
            comparison_results = pickle.load(f)
    print(
        f"[compare] comparison results for {len(comparison_results)} protein(s)",
        flush=True,
    )
    return detected_domains, comparison_results


def build_feature_matrix(
    ids: list[str],
    san_to_orig: dict[str, str],
    detected_domains: dict,
    comparison_results: dict,
    module_columns: list[str],
    module_to_type: dict[str, str],
) -> pd.DataFrame:
    """Build the per-protein 1 - TM-score domain-comparison matrix.

    For each protein and each reference module column, value = ``1 - best TM-score``
    over the protein's detected domains (matched by domain type) against that module,
    0.0 when no comparison exists. Replicates easy_predict.py's ``dom_feat`` fill +
    ``1 - dom_feat``, but over the full union of reference modules rather than a single
    fold's subset. The alpha1/alpha2 split mirrors easy_predict.py: the first detected
    alpha domain is alpha1, the second alpha2.

    Like easy_predict.py, a detected domain only fills reference-module columns of the
    SAME (alpha1/alpha2-resolved) domain type — even though Foldseek may return cross-type
    TM-scores for a given detected-domain PDB, the classifier indexes by
    ``domain_type_2_order_of_domain_modules[domain_type]`` only.
    """
    col_index = {m: i for i, m in enumerate(module_columns)}
    orig_to_san = {v: k for k, v in san_to_orig.items()}
    n_cols = len(module_columns)

    rows = []
    for orig in ids:
        feat = np.zeros(n_cols, dtype=np.float32)  # best TM-score so far (0 default)
        san = orig_to_san.get(orig)
        if san is not None and san in comparison_results and san in detected_domains:
            current_cmp = comparison_results[san]
            was_alpha_observed = False
            for det in detected_domains[san]:
                domain_type = det.domain
                detection_id = det.module_id
                if detection_id not in current_cmp:
                    continue
                known_to_tm = dict(current_cmp[detection_id])
                if domain_type == "alpha":
                    if not was_alpha_observed:
                        domain_type = "alpha1"
                        was_alpha_observed = True
                    else:
                        domain_type = "alpha2"
                # fill best TM-score per reference module of THIS domain type only
                for known_module_id, tm in known_to_tm.items():
                    ci = col_index.get(known_module_id)
                    if ci is None or module_to_type.get(known_module_id) != domain_type:
                        continue
                    if tm > feat[ci]:
                        feat[ci] = tm
        # easy_predict.py: dom_feat = 1 - dom_feat  (only over compared modules; the
        # zeros stay 1.0? — NO: easy_predict applies 1 - to the WHOLE per-classifier
        # vector, where uncompared entries are 0, giving 1.0). We faithfully reproduce
        # that: emit 1 - feat for every column.
        row = {"id": orig}
        inv = 1.0 - feat
        for m, ci in col_index.items():
            row[m] = float(inv[ci])
        rows.append(row)
    df = pd.DataFrame(rows, columns=["id"] + module_columns)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract EnzymeExplorer domain-comparison features (the production "
        "model's structure/function block) for a set of proteins with structures."
    )
    parser.add_argument(
        "--sequences_csv", required=True,
        help="CSV with the protein IDs (column --id_column). Defines the output row set.",
    )
    parser.add_argument("--id_column", default="Enzyme_marts_ID")
    parser.add_argument(
        "--structs_dir", required=True,
        help="Directory of <id>.pdb structures (ESMFold/AF), filename stem == id.",
    )
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--scratch_dir", default=None,
                        help="Working dir for staged structures + EE temp files.")
    parser.add_argument("--n_jobs", type=int, default=16)
    parser.add_argument(
        "--classifier_pkl", default="data/classifier_domain_and_plm_checkpoints.pkl",
        help="Production fold-classifier bundle (defines the reference-module columns).",
    )
    parser.add_argument(
        "--reuse_cached", action="store_true",
        help="Skip detection+comparison and reuse the cached pickles in the scratch dir "
        "(scratch/_temp/{detections,comparisons}.pkl); only rebuild the feature matrix.",
    )
    args = parser.parse_args()

    ids = pd.read_csv(args.sequences_csv)[args.id_column].astype(str).tolist()
    print(f"{len(ids)} input protein id(s) from {args.sequences_csv}", flush=True)

    with open(args.classifier_pkl, "rb") as f:
        fold_classifiers = pickle.load(f)
    module_columns, module_to_type = build_module_columns(fold_classifiers)
    print(f"{len(module_columns)} reference-module feature columns (union over folds)",
          flush=True)

    scratch = Path(args.scratch_dir or (Path(args.output_csv).parent / "_ee_domain_scratch"))
    scratch.mkdir(parents=True, exist_ok=True)

    san_structs_dir, san_to_orig, needed_df = stage_structures(args.structs_dir, ids, scratch)
    print(f"{len(san_to_orig)} structure(s) staged (sanitized) of {len(ids)} id(s)",
          flush=True)
    needed_csv = scratch / "needed_proteins.csv"
    needed_df.to_csv(needed_csv, index=False)

    detected_domains, comparison_results = run_detection_and_comparison(
        san_structs_dir, needed_csv, scratch, n_jobs=args.n_jobs,
        reuse_cached=args.reuse_cached,
    )

    df = build_feature_matrix(
        ids, san_to_orig, detected_domains, comparison_results,
        module_columns, module_to_type,
    )
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)

    n_with = sum(
        1 for orig in ids
        if {v: k for k, v in san_to_orig.items()}.get(orig) in comparison_results
    )
    print(
        f"Wrote {len(df)} rows x {len(module_columns)} domain features to "
        f"{args.output_csv}; {n_with} protein(s) had ≥1 domain comparison.",
        flush=True,
    )


if __name__ == "__main__":
    main()
