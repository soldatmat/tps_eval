"""Label-agnostic k-NN coarse-label transfer for tps_eval.

Predicts a coarse class for each generated TPS design by a distance-weighted vote
of its nearest MARTS-DB known-TPS neighbours, ensembled across the three similarity
spaces tps_eval already computes (max_sequence_identity, min_embedding_distance,
structural_identity). The label assignments are an INPUT (a ``reference_id,label``
CSV) — the logic is agnostic to what the labels mean (first-cyclization class,
size class, substrate, ...). Swap the label file to change the labeling.

Two entry points:

* ``transfer_labels``   — predict mode. Given the three per-design top-k CSVs (the
  ``<input>_<tool>_topk.csv`` files emitted by the tools' ``--top_k`` flag), a label
  file, and a calibration JSON, produce a CSV keyed by ``ID`` with the per-space and
  ensembled predicted label + calibrated confidence + nearest-neighbour similarity.

* ``calibrate``         — leave-one-out calibration on MARTS-DB. Given the three
  MARTS-DB **self** top-k CSVs (each query's neighbours exclude itself) and the label
  file, measure accuracy as a function of nearest-neighbour similarity, per space and
  ensembled, and write the calibration artifact (thresholds tau + similarity->P(correct)
  curve) consumed by predict mode.

Top-k CSV contract (all three tools, columns ``query_id,rank,neighbour_id,score``):
    max_sequence_identity   score = identity PERCENT in [0,100]   LARGER = closer
    min_embedding_distance  score = ESM-embedding L2 distance      SMALLER = closer
    structural_identity     score = foldseek TM-score in [0,1]     LARGER = closer

Foldseek splits multi-chain PDBs, so a structural ``neighbour_id`` may carry a
trailing ``_<chain>`` suffix (e.g. ``marts_E00123_A``) — it is stripped before
joining to the label file (see ``_strip_chain_suffix``).
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Spaces
# --------------------------------------------------------------------------- #
# The three similarity spaces tps_eval computes. ``larger_is_closer`` records the
# score direction; literature-prior tau seeds the calibration grid (and is the
# fallback tau when calibration has too few labeled neighbours above threshold):
#   ~40% sequence identity is the classic "twilight zone" floor for reliable
#   homology/function transfer (Rost 1999); TM-score ~0.5 is the fold-similarity
#   floor (Xu & Zhang 2010). Embedding distance has no universal prior, so its
#   tau is purely empirical (seeded at 0 = no floor).
SPACES: Tuple[str, ...] = ("sequence", "embedding", "structural")

SPACE_LARGER_IS_CLOSER: Dict[str, bool] = {
    "sequence": True,       # identity percent
    "embedding": False,     # L2 distance
    "structural": True,     # TM-score
}

# Literature-prior seeds for tau, expressed as a SIMILARITY in [0,1]
# (after the score->similarity conversion below).
SPACE_PRIOR_TAU: Dict[str, float] = {
    "sequence": 0.40,       # 40% identity -> similarity 0.40
    "embedding": 0.0,       # no prior floor; learned empirically
    "structural": 0.50,     # TM-score 0.5 fold floor
}


def score_to_similarity(space: str, score: float) -> float:
    """Map a raw per-space score to a similarity in [0, 1] (LARGER = closer).

    sequence:   identity percent / 100
    structural: TM-score (already in [0,1])
    embedding:  1 / (1 + L2 distance)  (distance 0 -> sim 1; distance->inf -> sim 0)
    """
    if score is None or (isinstance(score, float) and math.isnan(score)):
        return float("nan")
    if space == "sequence":
        return float(score) / 100.0
    if space == "structural":
        return float(score)
    if space == "embedding":
        return 1.0 / (1.0 + float(score))
    raise ValueError(f"unknown space {space!r}")


# --------------------------------------------------------------------------- #
# IO helpers
# --------------------------------------------------------------------------- #
def _strip_chain_suffix(neighbour_id: str, valid_ids: Optional[set] = None) -> str:
    """Strip a foldseek ``_<chain>`` suffix from a structural neighbour id.

    Foldseek splits multi-chain PDBs into per-chain targets named ``<stem>_<chain>``.
    The chain token is a short alphanumeric (usually a single letter/digit). We strip
    a trailing ``_<token>`` only when (a) the full id is NOT already a known label id
    and (b) the stripped id IS a known label id (when ``valid_ids`` is given), or
    (c) ``valid_ids`` is None and the trailing token looks like a chain (<=2 chars).
    This avoids clobbering legitimate ids that themselves contain underscores.
    """
    nid = str(neighbour_id)
    if valid_ids is not None:
        if nid in valid_ids:
            return nid
        if "_" in nid:
            stem = nid.rsplit("_", 1)[0]
            if stem in valid_ids:
                return stem
        return nid
    # No reference set to check against: strip a short trailing chain token.
    if "_" in nid:
        stem, tok = nid.rsplit("_", 1)
        if len(tok) <= 2 and tok.isalnum():
            return stem
    return nid


def load_label_map(label_file: str) -> Tuple[Dict[str, object], List[object]]:
    """Load the ``reference_id,label`` CSV into an id->label dict + sorted class list.

    Accepts a header with either ``reference_id,label`` or, more leniently, the first
    two columns as (id, label). Labels are kept as-is (str or int); the class list is
    the sorted set of distinct labels. ``unknown`` / NaN labels are dropped.
    """
    df = pd.read_csv(label_file)
    cols = list(df.columns)
    if "reference_id" in cols and "label" in cols:
        id_col, lab_col = "reference_id", "label"
    else:
        id_col, lab_col = cols[0], cols[1]
    label_map: Dict[str, object] = {}
    for rid, lab in zip(df[id_col].astype(str), df[lab_col]):
        if pd.isna(lab):
            continue
        label_map[rid] = lab
    classes = sorted({v for v in label_map.values()}, key=lambda x: str(x))
    return label_map, classes


def _topk_groups(topk_csv: str) -> Dict[str, List[Tuple[int, str, float]]]:
    """Read a top-k CSV into {query_id: [(rank, neighbour_id, score), ...]} (rank asc)."""
    df = pd.read_csv(topk_csv)
    df = df.sort_values(["query_id", "rank"], kind="stable")
    out: Dict[str, List[Tuple[int, str, float]]] = defaultdict(list)
    for qid, rank, nid, score in zip(
        df["query_id"].astype(str), df["rank"], df["neighbour_id"].astype(str), df["score"]
    ):
        out[qid].append((int(rank), nid, float(score)))
    return out


# --------------------------------------------------------------------------- #
# Per-space voting
# --------------------------------------------------------------------------- #
@dataclass
class SpaceVote:
    """Outcome of the distance-weighted vote in one space for one query."""

    predicted: Optional[object]      # argmax class, or None if abstained
    confidence: float                # raw (pre-calibration) confidence in [0,1]
    posterior: Dict[object, float]   # normalized class weights
    nn_similarity: float             # nearest-neighbour similarity in [0,1] (NaN if none)
    n_voters: int                    # neighbours above tau that contributed a label


def vote_space(
    neighbours: Sequence[Tuple[int, str, float]],
    space: str,
    label_map: Dict[str, object],
    classes: Sequence[object],
    tau: float,
    *,
    valid_ids: Optional[set] = None,
) -> SpaceVote:
    """Distance-weighted vote in one space.

    For each neighbour: convert score->similarity, strip chain suffix (structural),
    look up the label, and accumulate ``similarity`` as the vote weight — but only
    for neighbours whose similarity >= ``tau`` (others are ignored => abstain when
    none qualify). Posterior = normalized per-class weights; predicted = argmax.
    Confidence = winning_fraction * topk_agreement * nn_similarity, where
        winning_fraction = posterior[argmax]
        topk_agreement   = (# qualifying voters for the winning class) / (# qualifying voters)
        nn_similarity    = similarity of the single nearest qualifying neighbour.
    """
    weights: Dict[object, float] = defaultdict(float)
    counts: Dict[object, int] = defaultdict(int)
    nn_similarity = float("nan")
    best_sim = -1.0
    n_voters = 0
    for _rank, nid, score in neighbours:
        sim = score_to_similarity(space, score)
        if math.isnan(sim):
            continue
        key = _strip_chain_suffix(nid, valid_ids) if space == "structural" else nid
        lab = label_map.get(key)
        if lab is None:
            continue
        # Track nearest-neighbour similarity over labeled neighbours (any sim).
        if sim > best_sim:
            best_sim = sim
            nn_similarity = sim
        if sim < tau:
            continue  # below threshold -> does not vote
        weights[lab] += sim
        counts[lab] += 1
        n_voters += 1

    if n_voters == 0 or not weights:
        return SpaceVote(None, 0.0, {}, nn_similarity, 0)

    total = sum(weights.values())
    posterior = {c: weights.get(c, 0.0) / total for c in weights}
    predicted = max(posterior, key=posterior.get)
    winning_fraction = posterior[predicted]
    topk_agreement = counts[predicted] / n_voters
    confidence = winning_fraction * topk_agreement * max(0.0, nn_similarity)
    return SpaceVote(predicted, confidence, posterior, nn_similarity, n_voters)


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #
def _calibration_curve(
    similarities: np.ndarray, correct: np.ndarray, edges: Sequence[float]
) -> List[Dict[str, float]]:
    """Bin (nn_similarity -> P(correct)) accuracy curve over the given bin edges."""
    similarities = np.asarray(similarities, dtype=float)
    correct = np.asarray(correct, dtype=float)
    curve = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (similarities >= lo) & (similarities < hi)
        n = int(mask.sum())
        acc = float(correct[mask].mean()) if n > 0 else float("nan")
        curve.append({"lo": float(lo), "hi": float(hi), "n": n, "accuracy": acc})
    return curve


def _choose_tau(
    similarities: np.ndarray,
    correct: np.ndarray,
    prior_tau: float,
    target_accuracy: float = 0.5,
    min_support: int = 20,
) -> float:
    """Pick tau as the lowest nn_similarity at which empirical P(correct) >= target.

    Sweeps candidate thresholds; for each, computes accuracy over predictions whose
    nn_similarity >= candidate (requiring >= min_support such predictions). Returns
    the smallest candidate meeting ``target_accuracy``; falls back to ``prior_tau``
    if none qualifies (e.g. too little data). Never returns below ``prior_tau``.
    """
    similarities = np.asarray(similarities, dtype=float)
    correct = np.asarray(correct, dtype=float)
    finite = np.isfinite(similarities)
    similarities, correct = similarities[finite], correct[finite]
    if similarities.size == 0:
        return prior_tau
    candidates = np.linspace(0.0, 1.0, 101)
    best = None
    for c in candidates:
        mask = similarities >= c
        n = int(mask.sum())
        if n < min_support:
            continue
        acc = float(correct[mask].mean())
        if acc >= target_accuracy:
            best = float(c)
            break
    chosen = prior_tau if best is None else max(best, prior_tau)
    return float(chosen)


def calibrate(
    self_topk: Dict[str, str],
    label_file: str,
    *,
    top_k: Optional[int] = None,
    target_accuracy: float = 0.5,
    labeling: str = "labeling",
) -> dict:
    """Leave-one-out calibration on MARTS-DB.

    Args:
        self_topk: {space: path} to the MARTS-DB self top-k CSVs (self excluded).
        label_file: reference_id,label CSV.
        top_k: cap neighbours per query (default: use all present).
        target_accuracy: accuracy floor used to pick tau per space/ensemble.
        labeling: name recorded in the artifact (e.g. "first_cyclization").

    Returns the calibration dict (also the JSON artifact contents).
    """
    label_map, classes = load_label_map(label_file)
    valid_ids = set(label_map)

    # First pass: learn tau per space from raw nn_similarity-vs-correct, with a
    # provisional tau of the prior (so a neighbour only needs to be labeled to
    # contribute its nn_similarity, which is what tau is calibrated against).
    bin_edges = list(np.round(np.linspace(0.0, 1.0, 21), 4))
    calibration: dict = {
        "labeling": labeling,
        "classes": [str(c) for c in classes],
        "n_classes": len(classes),
        "spaces": {},
        "ensemble": {},
        "top_k": top_k,
        "target_accuracy": target_accuracy,
    }

    groups = {space: _topk_groups(path) for space, path in self_topk.items()}
    # Queries = labeled MARTS ids that appear in at least one space.
    query_ids = sorted(
        {q for g in groups.values() for q in g} & valid_ids
    )

    # --- learn per-space tau ---
    per_space_records: Dict[str, Dict[str, list]] = {}
    for space in self_topk:
        sims, correct = [], []
        g = groups[space]
        for qid in query_ids:
            true_lab = label_map[qid]
            nbrs = g.get(qid, [])
            if top_k is not None:
                nbrs = nbrs[:top_k]
            # tau=prior for the learning pass (just need labeled voters + nn sim)
            vote = vote_space(
                nbrs, space, label_map, classes, SPACE_PRIOR_TAU[space], valid_ids=valid_ids
            )
            if vote.predicted is None or math.isnan(vote.nn_similarity):
                continue
            sims.append(vote.nn_similarity)
            correct.append(1.0 if vote.predicted == true_lab else 0.0)
        sims_a, corr_a = np.array(sims), np.array(correct)
        tau = _choose_tau(sims_a, corr_a, SPACE_PRIOR_TAU[space], target_accuracy)
        per_space_records[space] = {"sims": sims_a, "correct": corr_a, "tau": tau}

    # --- per-space final accuracy + curve at the chosen tau ---
    for space in self_topk:
        rec = per_space_records[space]
        tau = rec["tau"]
        sims, correct = [], []
        g = groups[space]
        n_abstain = 0
        for qid in query_ids:
            true_lab = label_map[qid]
            nbrs = g.get(qid, [])
            if top_k is not None:
                nbrs = nbrs[:top_k]
            vote = vote_space(nbrs, space, label_map, classes, tau, valid_ids=valid_ids)
            if vote.predicted is None:
                n_abstain += 1
                continue
            sims.append(vote.nn_similarity)
            correct.append(1.0 if vote.predicted == true_lab else 0.0)
        sims_a, corr_a = np.array(sims), np.array(correct)
        calibration["spaces"][space] = {
            "tau": float(tau),
            "prior_tau": SPACE_PRIOR_TAU[space],
            "larger_is_closer": SPACE_LARGER_IS_CLOSER[space],
            "n_predicted": int(len(corr_a)),
            "n_abstained": int(n_abstain),
            "n_queries": int(len(query_ids)),
            "accuracy": float(corr_a.mean()) if corr_a.size else float("nan"),
            "calibration_curve": _calibration_curve(sims_a, corr_a, bin_edges),
        }

    # --- ensemble: average calibrated posteriors across spaces ---
    taus = {space: per_space_records[space]["tau"] for space in self_topk}
    ens_sims, ens_correct = [], []
    n_abstain = 0
    for qid in query_ids:
        true_lab = label_map[qid]
        agg: Dict[object, float] = defaultdict(float)
        n_contrib = 0
        max_nn = float("nan")
        for space in self_topk:
            nbrs = groups[space].get(qid, [])
            if top_k is not None:
                nbrs = nbrs[:top_k]
            vote = vote_space(nbrs, space, label_map, classes, taus[space], valid_ids=valid_ids)
            if vote.predicted is None:
                continue
            n_contrib += 1
            for c, p in vote.posterior.items():
                agg[c] += p
            if math.isnan(max_nn) or (not math.isnan(vote.nn_similarity) and vote.nn_similarity > max_nn):
                max_nn = vote.nn_similarity
        if n_contrib == 0 or not agg:
            n_abstain += 1
            continue
        pred = max(agg, key=agg.get)
        ens_sims.append(max_nn)
        ens_correct.append(1.0 if pred == true_lab else 0.0)
    ens_sims_a, ens_corr_a = np.array(ens_sims), np.array(ens_correct)
    calibration["ensemble"] = {
        "n_predicted": int(len(ens_corr_a)),
        "n_abstained": int(n_abstain),
        "n_queries": int(len(query_ids)),
        "accuracy": float(ens_corr_a.mean()) if ens_corr_a.size else float("nan"),
        "calibration_curve": _calibration_curve(ens_sims_a, ens_corr_a, bin_edges),
    }
    return calibration


def _interp_calibrated_confidence(nn_similarity: float, curve: Sequence[dict]) -> float:
    """Map an nn_similarity to calibrated P(correct) via the binned curve.

    Uses the bin containing nn_similarity; if that bin has no support, walk down to
    the nearest lower populated bin. Returns NaN if no populated bin at/below.
    """
    if math.isnan(nn_similarity):
        return float("nan")
    populated = [b for b in curve if b["n"] > 0 and not math.isnan(b["accuracy"])]
    if not populated:
        return float("nan")
    # exact bin
    for b in curve:
        if b["lo"] <= nn_similarity < b["hi"] and b["n"] > 0 and not math.isnan(b["accuracy"]):
            return float(b["accuracy"])
    # nearest lower populated bin
    lower = [b for b in populated if b["lo"] <= nn_similarity]
    if lower:
        return float(max(lower, key=lambda b: b["lo"])["accuracy"])
    # below all -> lowest populated
    return float(min(populated, key=lambda b: b["lo"])["accuracy"])


# --------------------------------------------------------------------------- #
# Predict
# --------------------------------------------------------------------------- #
ABSTAIN_LABEL = "unknown"


def transfer_labels(
    design_topk: Dict[str, str],
    label_file: str,
    calibration: dict,
    *,
    top_k: Optional[int] = None,
) -> pd.DataFrame:
    """Predict a coarse label per design and return a DataFrame keyed by ``ID``.

    Args:
        design_topk: {space: path} to the per-design top-k CSVs. Missing spaces are
            simply skipped (the design abstains in that space).
        label_file: reference_id,label CSV (the same labeling the calibration used).
        calibration: loaded calibration dict (provides tau + curves per space/ensemble).
        top_k: cap neighbours per query (default: all present).

    Output columns (keyed by ``ID``):
        predicted_label, confidence,
        predicted_label_<space>, conf_<space>, nn_similarity_<space>  (each space).
    Designs below tau in ALL spaces abstain: predicted_label = "unknown", confidence
    low (NaN-derived -> 0). Novel designs SHOULD land here.
    """
    label_map, classes = load_label_map(label_file)
    valid_ids = set(label_map)
    spaces = [s for s in SPACES if s in design_topk]

    groups = {space: _topk_groups(design_topk[space]) for space in spaces}
    all_ids = sorted({q for g in groups.values() for q in g})
    taus = {s: calibration["spaces"][s]["tau"] for s in spaces}
    curves = {s: calibration["spaces"][s]["calibration_curve"] for s in spaces}
    ens_curve = calibration.get("ensemble", {}).get("calibration_curve", [])

    rows = []
    for qid in all_ids:
        row: Dict[str, object] = {"ID": qid}
        agg: Dict[object, float] = defaultdict(float)
        n_contrib = 0
        max_nn = float("nan")
        for space in spaces:
            nbrs = groups[space].get(qid, [])
            if top_k is not None:
                nbrs = nbrs[:top_k]
            vote = vote_space(nbrs, space, label_map, classes, taus[space], valid_ids=valid_ids)
            row[f"predicted_label_{space}"] = (
                str(vote.predicted) if vote.predicted is not None else ABSTAIN_LABEL
            )
            row[f"nn_similarity_{space}"] = vote.nn_similarity
            # calibrated per-space confidence from the space's curve
            if vote.predicted is None:
                row[f"conf_{space}"] = 0.0
            else:
                cal = _interp_calibrated_confidence(vote.nn_similarity, curves[space])
                row[f"conf_{space}"] = cal if not math.isnan(cal) else vote.confidence
                n_contrib += 1
                for c, p in vote.posterior.items():
                    agg[c] += p
                if math.isnan(max_nn) or (
                    not math.isnan(vote.nn_similarity) and vote.nn_similarity > max_nn
                ):
                    max_nn = vote.nn_similarity
        if n_contrib == 0 or not agg:
            row["predicted_label"] = ABSTAIN_LABEL
            row["confidence"] = 0.0
        else:
            pred = max(agg, key=agg.get)
            row["predicted_label"] = str(pred)
            ens_conf = _interp_calibrated_confidence(max_nn, ens_curve)
            if math.isnan(ens_conf):
                # fall back to mean of per-space calibrated confidences
                confs = [row[f"conf_{s}"] for s in spaces if row.get(f"predicted_label_{s}") != ABSTAIN_LABEL]
                ens_conf = float(np.mean(confs)) if confs else 0.0
            row["confidence"] = float(ens_conf)
        rows.append(row)

    # Stable column order.
    cols = ["ID", "predicted_label", "confidence"]
    for space in spaces:
        cols += [f"predicted_label_{space}", f"conf_{space}", f"nn_similarity_{space}"]
    return pd.DataFrame(rows, columns=cols)


def save_calibration(calibration: dict, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(calibration, fh, indent=2)


def load_calibration(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)
