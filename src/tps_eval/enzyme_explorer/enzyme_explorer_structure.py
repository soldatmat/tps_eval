"""Structure-based TPS classification per design, via EnzymeExplorer.

Reshapes the TWO CSVs written by EnzymeExplorer's ``predict_with_structures``
console script into ONE CSV keyed by ``ID``, one row per input design, so the
structure-mode isTPS score can be consumed like any other tps_eval metric.

Why this exists — and why "no row" is the whole point
-----------------------------------------------------
``predict_with_structures`` runs two classifiers and writes one CSV per
classifier, because each has its own calibrator and their scores are therefore
NOT comparable:

  * ``predictions_plm_domains.csv``      — the PLM_Domains classifier. Scores
    only those designs whose structure yielded at least one detected TPS
    structural domain (alpha/beta/gamma/ids/delta/epsilon/zeta).
  * ``predictions_plm_only_fallback.csv`` — the PLM-only classifier. Scores the
    REST: designs for which the domain detector returned nothing. EE routes them
    here in ``ensemble.predict_with_plm_and_domains`` (a protein absent from
    ``structural_features_ids`` lands in ``fallback_ids``).

So a design with NO detected TPS domain simply never appears in the domains CSV.
That absence is a MEANINGFUL REJECTION — "this fold has no TPS-like domain
architecture" — not a missing datum. Accordingly:

  * ``ee_struct_status`` records which of the two happened
    (``scored`` / ``no_domains`` / ``missing``);
  * ``isTPS_struct_raw`` / ``isTPS_struct_cal`` are NaN for a ``no_domains``
    design, so any downstream ``>= threshold`` test evaluates False and the
    design FAILS the structure-mode filter;
  * the PLM-only fallback score is carried in a SEPARATELY NAMED column
    (``isTPS_fallback_raw`` / ``_cal``) for diagnostics ONLY. Do NOT substitute
    it for the structure score — it is a different classifier with a different
    calibration and it is exactly the score that gets fooled by degenerate
    sequences.

⚠ Run the prefilter OFF for any reported number. ``predict_with_structures``
takes ``--prefilter-pdbs-by-foldseek``, an opt-in speedup whose own docstring
admits a "small recall loss" — for a GATE, recall loss means falsely rejecting
real designs. EnzymeExplorer's ``revision`` branch HARDCODES it on
(``prediction/domains.py``) with no CLI override; the ``main`` branch defaults it
off. tps_eval therefore points ``ee_struct`` at the ``main`` checkout via
``ENZYME_EXPLORER_STRUCT_{PATH,ENV}`` — see the comment in ``paths.sh``. This
matters concretely: in the one pre-2026-08-26 run (job 1488255) only 1 of 56
structures cleared the prefilter, so domain detection never ran on the other 55
and its "0/56 domain hits" headline could not be separated from a prefilter
artifact.

The full design universe is enumerated from the INPUT sequences, not from either
output CSV, and left-joined — so a design EE dropped entirely (e.g. length) gets
a row with ``ee_struct_status="missing"`` rather than vanishing from the metric.

Column names here are deliberately INDEPENDENT of EnzymeExplorer's own schema,
which differs between checkouts: the ``main`` branch emits ``<class>_raw`` /
``<class>_p`` while the older ``revision`` branch emits ``<class>_score`` /
``<class>_p_calibrated``. Both are accepted on read.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# EnzymeExplorer's fixed class order (see its prediction/calibration.py). "TPS"
# is the is-it-a-terpene-synthase call; the rest are per-substrate classes.
EE_CLASSES: List[str] = [
    "TPS", "GPP", "FPP", "GGPP", "GFPP", "CPP", "EDSQ", "2xFPP", "2xGGPP", "IDS",
]

# Per-design status values, in the order they are documented above.
STATUS_SCORED = "scored"        # >=1 TPS domain detected -> PLM_Domains score
STATUS_NO_DOMAINS = "no_domains"  # zero domains detected -> PLM-only fallback
STATUS_MISSING = "missing"      # absent from BOTH CSVs (EE dropped it)

DOMAINS_CSV = "predictions_plm_domains.csv"
FALLBACK_CSV = "predictions_plm_only_fallback.csv"

_RAW_SUFFIXES = ("_raw", "_score")          # main branch, revision branch
_CAL_SUFFIXES = ("_p", "_p_calibrated")     # main branch, revision branch

COLUMNS: List[str] = (
    ["ID", "ee_struct_status", "ee_struct_domain_hit"]
    + [f"{c}_struct_raw" for c in EE_CLASSES]
    + [f"{c}_struct_cal" for c in EE_CLASSES]
    + ["isTPS_fallback_raw", "isTPS_fallback_cal"]
)


def _pick_column(df: pd.DataFrame, klass: str, suffixes: Tuple[str, ...]) -> Optional[str]:
    """First column named ``<klass><suffix>`` present in `df`, else None.

    Accepts either EnzymeExplorer schema (see the module docstring). ``_p`` is
    tried before ``_p_calibrated`` and ``_raw`` before ``_score``, matching the
    branch we run; a checkout emitting only the other naming still resolves.
    """
    for suffix in suffixes:
        col = f"{klass}{suffix}"
        if col in df.columns:
            return col
    return None


def _read_predictions(path: Path) -> pd.DataFrame:
    """Read one of EE's prediction CSVs, normalised to ``ID`` + per-class floats.

    An EE run in which every design took the SAME branch leaves the other CSV
    degenerate, and it does so in two different shapes — both legitimate
    outcomes, neither a failure:

      * **header only, zero rows** — the frame had columns but no records
        (observed as ``id,sequence`` when the empty frame kept only the
        passthrough columns);
      * **a 1-byte file containing just a newline** — the frame had no columns
        at all, so ``to_csv`` wrote nothing. ``pd.read_csv`` raises
        ``EmptyDataError`` on this, which is why it must be caught rather than
        allowed to abort the arm. This is the COMMON case, not an edge case: it
        happens whenever every design gets a domain hit (verified on 300 natural
        MARTS-DB TPS, where 300/300 scored and the fallback CSV was 1 byte).
    """
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["ID"])
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["ID"])
    if df.empty:
        return pd.DataFrame(columns=["ID"])
    id_col = "id" if "id" in df.columns else ("ID" if "ID" in df.columns else None)
    if id_col is None:
        raise ValueError(f"{path}: no 'id'/'ID' column (columns: {list(df.columns)})")
    out = pd.DataFrame({"ID": df[id_col].astype(str)})
    for klass in EE_CLASSES:
        raw_col = _pick_column(df, klass, _RAW_SUFFIXES)
        cal_col = _pick_column(df, klass, _CAL_SUFFIXES)
        out[f"{klass}_raw"] = (
            pd.to_numeric(df[raw_col], errors="coerce") if raw_col else float("nan")
        )
        out[f"{klass}_cal"] = (
            pd.to_numeric(df[cal_col], errors="coerce") if cal_col else float("nan")
        )
    return out


def design_ids_from_fasta(fasta_path: str) -> List[str]:
    """Every design ID in a FASTA, in file order (ID == first token of the header)."""
    ids: List[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def design_ids_from_csv(csv_path: str, id_column: str = "id") -> List[str]:
    """Every design ID in a sequences CSV, in file order."""
    df = pd.read_csv(csv_path)
    col = id_column if id_column in df.columns else ("ID" if "ID" in df.columns else None)
    if col is None:
        raise ValueError(f"{csv_path}: no '{id_column}'/'ID' column")
    return [str(x) for x in df[col]]


def design_ids_from_structs_dir(structs_dir: str) -> List[str]:
    """Every structure's ID (== ``.pdb`` filename stem) in `structs_dir`.

    Only a fallback ID universe: EE's detector globs ``*.pdb``, so this misses
    designs that never folded — and those SHOULD appear in the metric as
    failures. Prefer the input sequences when available.
    """
    return sorted({p.stem for p in Path(structs_dir).glob("*.pdb")})


def build_structure_predictions(
    output_dir: str,
    design_ids: List[str],
) -> pd.DataFrame:
    """Reshape EE's two prediction CSVs in `output_dir` into one per-design frame.

    `design_ids` is the authoritative design universe; every one of them gets
    exactly one row. See the module docstring for the status semantics.
    """
    out_dir = Path(output_dir)
    domains = _read_predictions(out_dir / DOMAINS_CSV)
    fallback = _read_predictions(out_dir / FALLBACK_CSV)

    overlap = set(domains["ID"]) & set(fallback["ID"])
    if overlap:
        # The two passes partition the input by construction; an overlap means we
        # are reading mismatched CSVs (e.g. two runs written into one dir), and
        # silently keeping one score would misreport the gate.
        raise ValueError(
            f"{output_dir}: {len(overlap)} design(s) appear in BOTH "
            f"{DOMAINS_CSV} and {FALLBACK_CSV} (e.g. {sorted(overlap)[:3]}); "
            "the two passes must partition the input — is this a stale/mixed "
            "output directory?"
        )

    df = pd.DataFrame({"ID": [str(i) for i in design_ids]})
    if df["ID"].duplicated().any():
        dupes = df.loc[df["ID"].duplicated(), "ID"].unique()[:3]
        raise ValueError(f"duplicate design IDs in input (e.g. {list(dupes)})")

    scored = set(domains["ID"])
    fell_back = set(fallback["ID"])
    df["ee_struct_status"] = [
        STATUS_SCORED if i in scored else (STATUS_NO_DOMAINS if i in fell_back else STATUS_MISSING)
        for i in df["ID"]
    ]
    df["ee_struct_domain_hit"] = df["ee_struct_status"] == STATUS_SCORED

    dom_idx = domains.set_index("ID") if not domains.empty else None
    fb_idx = fallback.set_index("ID") if not fallback.empty else None

    for klass in EE_CLASSES:
        for flavour, ee_suffix in (("raw", "_raw"), ("cal", "_cal")):
            col = f"{klass}_struct_{flavour}"
            src = f"{klass}{ee_suffix}"
            df[col] = (
                df["ID"].map(dom_idx[src]) if dom_idx is not None and src in dom_idx.columns
                else float("nan")
            )

    # Diagnostics only -- a DIFFERENT classifier with its own calibration. Never
    # a substitute for the structure score (see the module docstring).
    for flavour, ee_suffix in (("raw", "_raw"), ("cal", "_cal")):
        col = f"isTPS_fallback_{flavour}"
        src = f"TPS{ee_suffix}"
        df[col] = (
            df["ID"].map(fb_idx[src]) if fb_idx is not None and src in fb_idx.columns
            else float("nan")
        )

    # Canonical aliases for the TPS class -- the gate's actual input.
    df["isTPS_struct_raw"] = df["TPS_struct_raw"]
    df["isTPS_struct_cal"] = df["TPS_struct_cal"]

    ordered = ["ID", "ee_struct_status", "ee_struct_domain_hit",
               "isTPS_struct_raw", "isTPS_struct_cal"]
    ordered += [c for c in COLUMNS if c not in ordered]
    return df[ordered].sort_values("ID").reset_index(drop=True)


def load_id_map(id_map_csv: str) -> Dict[str, str]:
    """Read a ``sanitized_id,original_id`` CSV into a mapping.

    Needed because EnzymeExplorer's ``get_pdb_files`` REJECTS any structure whose
    filename stem contains a character outside ``[a-zA-Z0-9_]`` -- it raises on the
    first offender, aborting the whole directory. Some of our arms have dots or
    hyphens baked into every design ID (e.g. a sweep cell named ``foldcurve_s0.5``),
    so those runs go through a sanitized symlink farm and the original IDs are
    restored here, rather than being silently dropped from the metric.
    """
    df = pd.read_csv(id_map_csv)
    missing = {"sanitized_id", "original_id"} - set(df.columns)
    if missing:
        raise ValueError(f"{id_map_csv}: missing column(s) {sorted(missing)}")
    mapping = dict(zip(df["sanitized_id"].astype(str), df["original_id"].astype(str)))
    if len(mapping) != len(df):
        raise ValueError(f"{id_map_csv}: duplicate sanitized_id values")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f"{id_map_csv}: sanitization collided -- two IDs map to one")
    return mapping


def default_save_path(structs_dir: str) -> str:
    """``<structs_dir>_enzyme_explorer_structure.csv``, a sibling of the dir.

    Keyed off the structures directory (not the FASTA) to match every other
    structure-branch metric in the pipeline.
    """
    d = structs_dir.rstrip(os.sep)
    return os.path.join(
        os.path.dirname(d), os.path.basename(d) + "_enzyme_explorer_structure.csv"
    )


def reshape_to_csv(
    output_dir: str,
    *,
    structs_dir: str,
    fasta_path: Optional[str] = None,
    sequences_csv_path: Optional[str] = None,
    csv_id_column: str = "id",
    save_path: Optional[str] = None,
    id_map_csv: Optional[str] = None,
) -> str:
    """Write the per-design structure-mode CSV; returns the path written.

    The design universe comes from `fasta_path` / `sequences_csv_path` when
    given (the true set of designs, including any that never folded), else from
    the ``.pdb`` stems in `structs_dir`. `id_map_csv` restores original IDs after a
    sanitized run (see `load_id_map`).
    """
    if fasta_path:
        design_ids = design_ids_from_fasta(fasta_path)
    elif sequences_csv_path:
        design_ids = design_ids_from_csv(sequences_csv_path, csv_id_column)
    else:
        design_ids = design_ids_from_structs_dir(structs_dir)

    df = build_structure_predictions(output_dir, design_ids)
    if id_map_csv:
        # The run used sanitized IDs; restore the originals so the CSV joins against
        # every other per-design metric. Refuse rather than half-rename: a partial map
        # would silently produce a table keyed by two different ID conventions.
        mapping = load_id_map(id_map_csv)
        unmapped = sorted(set(df["ID"]) - set(mapping))
        if unmapped:
            raise ValueError(
                f"{id_map_csv} does not cover {len(unmapped)} ID(s) present in the "
                f"predictions (e.g. {unmapped[:3]}); refusing to write a table with "
                "mixed ID conventions")
        df["ID"] = df["ID"].map(mapping)
        df = df.sort_values("ID").reset_index(drop=True)
    path = save_path or default_save_path(structs_dir)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def summarize(df: pd.DataFrame, threshold: float = 0.5) -> Dict[str, int]:
    """Counts for a run log: how the design universe split across the branches."""
    n = len(df)
    hit = int(df["ee_struct_domain_hit"].sum())
    return {
        "n_designs": n,
        "n_domain_hit": hit,
        "n_no_domains": int((df["ee_struct_status"] == STATUS_NO_DOMAINS).sum()),
        "n_missing": int((df["ee_struct_status"] == STATUS_MISSING).sum()),
        "n_isTPS_struct_pass": int((df["isTPS_struct_raw"] >= threshold).sum()),
    }
