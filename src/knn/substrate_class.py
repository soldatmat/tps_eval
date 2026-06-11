"""Substrate-class combiner for tps_eval — predict each design's prenyl-diphosphate
substrate class (GPP/C10 mono, FPP/C15 sesqui, GGPP/C20 di, ...) by FUSING three
existing tps_eval signals. It is a THIN combiner: it computes nothing new about the
structures or sequences, it only reconciles outputs other tools already produced.

Signals fused (each is optional; the combiner uses whatever is present):

1. **k-NN substrate vote** (primary) — `knn_label_transfer.transfer_labels` run with the
   SUBSTRATE label file (`substrate_labels.csv`, MARTS `Type` -> GPP/FPP/GGPP/...) and the
   substrate calibration (`knn_calibration_substrate.json`). A distance-weighted vote of
   the design's nearest MARTS-DB TPS neighbours across the sequence/embedding/structural
   spaces -> `predicted_substrate` + calibrated `confidence`. This is the call of record.

2. **Pocket volume** (`pocket_descriptors.catalytic_pocket_volume`) — the active-site
   cavity is the molecular ruler that selects substrate chain length, so volume tracks
   substrate size monotonically (mono < sesqui < di < ...). We map the raw fpocket volume
   to a coarse SIZE band and check whether that band's size-rank is CONSISTENT with the
   k-NN call (`substrate_agreement`). fpocket absolute volumes are noisy/instrument-
   dependent, so the bands are deliberately coarse and the agreement is a "within one
   size class" sanity check, not a hard vote. (See POCKET_VOLUME_BANDS for the thresholds
   and their provenance.)

3. **EnzymeExplorer sequence-only** (`<input>_enzyme_explorer_sequence_only.csv`) — EE
   emits a per-substrate probability for each candidate prenyl-PP (`FPP_score`, `GPP_score`,
   `GGPP_score`, `GFPP_score`, `EDSQ_score`, ...; columns chosen so they share the label
   vocabulary that `make_substrate_labels.py` assigns). We take EE's argmax substrate as an
   independent classifier (`ee_substrate` + `ee_score`) and whether it agrees with k-NN.

Output `<input>_substrate_class.csv` keyed by `ID`:
    predicted_substrate    final substrate class (the k-NN call; "unknown" if k-NN abstains
                           AND no fallback signal resolves it)
    confidence             k-NN calibrated ensemble confidence in [0,1]
    knn_substrate          raw k-NN call (may be "unknown")
    knn_confidence         k-NN confidence
    pocket_volume          catalytic_pocket_volume (A^3), NaN if unavailable
    pocket_volume_band     coarse substrate-size band from pocket volume (or "" if NaN)
    substrate_agreement    do the k-NN call and the pocket-volume band agree (within one
                           size class)? True/False/NA (NA when either signal is missing)
    ee_substrate           EnzymeExplorer argmax substrate (or "" if EE not provided)
    ee_score               EE score of that substrate
    ee_agreement           does EE's argmax match the k-NN call? True/False/NA
    n_signals_agree        how many of {pocket-band, EE} corroborate the k-NN call

When k-NN abstains ("unknown") but EE is confident, `predicted_substrate` falls back to the
EE argmax (flagged via `predicted_substrate_source`).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from knn_label_transfer import (  # noqa: E402
    ABSTAIN_LABEL,
    load_calibration,
    transfer_labels,
)

# --------------------------------------------------------------------------- #
# Substrate size ordering (carbon count). Used to (a) rank pocket-volume bands
# and (b) measure "within one size class" agreement.
# --------------------------------------------------------------------------- #
SUBSTRATE_CARBONS: Dict[str, int] = {
    "DMAPP": 5,
    "GPP": 10,
    "FPP": 15,
    "GGPP": 20,
    "GFPP": 25,
    "EDSQ": 30,
    "C35": 35,
    "2xGGPP": 40,
    # IDS = prenyltransferase / chain elongation: no single substrate size -> excluded
    # from the size-ordered agreement checks (treated as "no size rank").
}

# Size-ordered list of the substrate classes that HAVE a meaningful chain length.
_SIZE_ORDER: List[str] = sorted(SUBSTRATE_CARBONS, key=lambda s: SUBSTRATE_CARBONS[s])
_SIZE_RANK: Dict[str, int] = {s: i for i, s in enumerate(_SIZE_ORDER)}

# --------------------------------------------------------------------------- #
# Pocket-volume -> substrate-size band.
# --------------------------------------------------------------------------- #
# The class-I TPS active-site cavity is the "molecular ruler" that gates substrate chain
# length, so the catalytic-pocket volume scales with substrate size (Christianson 2017,
# Chem. Rev.; Lesburg et al. 1997). fpocket absolute volumes are noisy and engine-
# dependent, so these are COARSE, conservative thresholds (A^3) chosen to separate the
# common mono/sesqui/di regime; they only feed a "within one size class" consistency
# check, never a hard vote. Each (upper_bound_exclusive, band_label): the band is the
# first whose volume < upper bound. Provenance: heuristic, to be replaced by a data-driven
# band once MARTS pocket-volume reference stats exist (none aggregated yet).
POCKET_VOLUME_BANDS: Tuple[Tuple[float, str], ...] = (
    (400.0, "GPP"),       # small cavity -> C10 monoterpene
    (900.0, "FPP"),       # mid cavity   -> C15 sesquiterpene
    (1600.0, "GGPP"),     # larger       -> C20 diterpene
    (float("inf"), "GFPP"),  # largest    -> C25+ (sester / larger)
)


def pocket_volume_to_band(volume: float) -> Optional[str]:
    """Map a raw catalytic-pocket volume (A^3) to a coarse substrate-size band, or None
    if the volume is NaN/unavailable."""
    if volume is None or (isinstance(volume, float) and np.isnan(volume)):
        return None
    for upper, band in POCKET_VOLUME_BANDS:
        if float(volume) < upper:
            return band
    return POCKET_VOLUME_BANDS[-1][1]


def _within_one_size_class(a: str, b: str) -> Optional[bool]:
    """Do substrate classes `a` and `b` agree within one size class? None if either has
    no size rank (e.g. IDS) — agreement is undefined for non-size-ranked classes."""
    if a not in _SIZE_RANK or b not in _SIZE_RANK:
        return None
    return abs(_SIZE_RANK[a] - _SIZE_RANK[b]) <= 1


# --------------------------------------------------------------------------- #
# EnzymeExplorer sequence-only signal
# --------------------------------------------------------------------------- #
# EE seq-only emits "<SMILES> (<name>)" columns in the structured output, but with the
# `predict_sequences_only` console script the columns are short codes:
#   2xFPP_score 2xGGPP_score CPP_score EDSQ_score FPP_score GFPP_score GGPP_score
#   GPP_score IDS_score TPS_score  (+ *_p_calibrated, isTPS, ID).
# We read the *_score columns, drop the non-substrate TPS/isTPS gate, fold CPP (copalyl-
# PP, a C20 diterpene precursor) into GGPP and 2xFPP into EDSQ so EE's vocabulary matches
# the substrate-label vocabulary, and take the argmax as EE's substrate call.
_EE_SCORE_SUFFIX = "_score"
_EE_NON_SUBSTRATE = {"TPS", "isTPS"}
# Fold EE-specific codes onto the shared substrate-label vocabulary.
_EE_CLASS_FOLD = {
    "CPP": "GGPP",    # copalyl-PP -> C20 diterpene
    "2xFPP": "EDSQ",  # 2xFPP -> squalene / epoxysqualene (C30)
}


def _ee_substrate_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Map substrate-class -> EE score column for the per-substrate *_score columns,
    skipping the TPS/isTPS gate columns and folding EE codes onto the shared vocabulary."""
    out: Dict[str, str] = {}
    for col in df.columns:
        c = str(col).strip()
        if not c.endswith(_EE_SCORE_SUFFIX):
            continue
        code = c[: -len(_EE_SCORE_SUFFIX)]
        if code in _EE_NON_SUBSTRATE:
            continue
        out[_EE_CLASS_FOLD.get(code, code)] = col
    return out


def load_ee_substrate(ee_csv: str) -> Dict[str, Tuple[str, float]]:
    """Read EE seq-only output into {ID: (argmax_substrate, score)}.

    Folded classes (CPP->GGPP, 2xFPP->EDSQ) are merged by taking the max score among the
    EE codes mapping to the same substrate class. Returns {} if the file has no usable
    per-substrate score columns.
    """
    df = pd.read_csv(ee_csv)
    df.columns = [str(c).strip() for c in df.columns]
    id_col = "ID" if "ID" in df.columns else ("id" if "id" in df.columns else df.columns[-1])
    cols = _ee_substrate_columns(df)
    if not cols:
        return {}
    out: Dict[str, Tuple[str, float]] = {}
    for _, row in df.iterrows():
        scores: Dict[str, float] = {}
        for sub, col in cols.items():
            try:
                v = float(row[col])
            except (TypeError, ValueError):
                continue
            if np.isnan(v):
                continue
            scores[sub] = max(scores.get(sub, float("-inf")), v)
        if not scores:
            continue
        argmax = max(scores, key=scores.get)
        out[str(row[id_col]).strip()] = (argmax, float(scores[argmax]))
    return out


def load_pocket_volume(pocket_csv: str) -> Dict[str, float]:
    """Read {ID: catalytic_pocket_volume} from a pocket_descriptors CSV."""
    df = pd.read_csv(pocket_csv)
    if "ID" not in df.columns or "catalytic_pocket_volume" not in df.columns:
        raise ValueError(
            f"{pocket_csv} must have columns ID + catalytic_pocket_volume "
            f"(got {list(df.columns)})."
        )
    out: Dict[str, float] = {}
    for _, row in df.iterrows():
        out[str(row["ID"]).strip()] = float(row["catalytic_pocket_volume"])
    return out


# --------------------------------------------------------------------------- #
# Combiner
# --------------------------------------------------------------------------- #
OUTPUT_COLUMNS = [
    "ID",
    "predicted_substrate",
    "predicted_substrate_source",
    "confidence",
    "knn_substrate",
    "knn_confidence",
    "pocket_volume",
    "pocket_volume_band",
    "substrate_agreement",
    "ee_substrate",
    "ee_score",
    "ee_agreement",
    "n_signals_agree",
]


def combine_substrate_class(
    design_topk: Dict[str, str],
    label_file: str,
    calibration: dict,
    *,
    pocket_csv: Optional[str] = None,
    ee_csv: Optional[str] = None,
    top_k: Optional[int] = None,
) -> pd.DataFrame:
    """Fuse the k-NN substrate vote with the pocket-volume band and the EE substrate signal.

    Args:
        design_topk: {space: path} per-design top-k CSVs (the k-NN inputs).
        label_file: SUBSTRATE reference_id,label CSV (substrate_labels.csv).
        calibration: loaded substrate calibration dict (knn_calibration_substrate.json).
        pocket_csv: optional <structs_dir>_pocket_descriptors.csv (catalytic_pocket_volume).
        ee_csv: optional <input>_enzyme_explorer_sequence_only.csv.
        top_k: cap neighbours per query.

    Returns a DataFrame keyed by ID (OUTPUT_COLUMNS).
    """
    knn = transfer_labels(design_topk, label_file, calibration, top_k=top_k)
    knn = knn.set_index("ID")

    pocket = load_pocket_volume(pocket_csv) if pocket_csv else {}
    ee = load_ee_substrate(ee_csv) if ee_csv else {}

    # All IDs seen in any signal.
    ids = set(knn.index.astype(str)) | set(pocket) | set(ee)

    rows: List[Dict[str, object]] = []
    for qid in sorted(ids):
        knn_sub = ABSTAIN_LABEL
        knn_conf = 0.0
        if qid in knn.index:
            knn_sub = str(knn.loc[qid, "predicted_label"])
            knn_conf = float(knn.loc[qid, "confidence"])

        # pocket-volume band
        vol = pocket.get(qid, float("nan"))
        band = pocket_volume_to_band(vol)

        # EE argmax
        ee_sub, ee_score = ("", float("nan"))
        if qid in ee:
            ee_sub, ee_score = ee[qid]

        # agreements (vs the k-NN call, when k-NN didn't abstain)
        have_knn = knn_sub != ABSTAIN_LABEL
        if have_knn and band is not None:
            agree = _within_one_size_class(knn_sub, band)
            substrate_agreement = bool(agree) if agree is not None else "NA"
        else:
            substrate_agreement = "NA"

        if have_knn and ee_sub:
            ee_agreement = bool(ee_sub == knn_sub)
        else:
            ee_agreement = "NA"

        n_agree = sum(1 for a in (substrate_agreement, ee_agreement) if a is True)

        # final call: k-NN of record; fall back to EE argmax when k-NN abstains
        if have_knn:
            predicted = knn_sub
            source = "knn"
            confidence = knn_conf
        elif ee_sub:
            predicted = ee_sub
            source = "enzyme_explorer"
            confidence = float(ee_score) if not np.isnan(ee_score) else 0.0
        else:
            predicted = ABSTAIN_LABEL
            source = "none"
            confidence = 0.0

        rows.append(
            {
                "ID": qid,
                "predicted_substrate": predicted,
                "predicted_substrate_source": source,
                "confidence": confidence,
                "knn_substrate": knn_sub,
                "knn_confidence": knn_conf,
                "pocket_volume": vol,
                "pocket_volume_band": band if band is not None else "",
                "substrate_agreement": substrate_agreement,
                "ee_substrate": ee_sub,
                "ee_score": ee_score,
                "ee_agreement": ee_agreement,
                "n_signals_agree": n_agree,
            }
        )

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
