from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Sequence, Tuple

import os
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers


SUBSTITUTION_MATRIX = substitution_matrices.load("BLOSUM62")

# Module-level aligner: global Needleman-Wunsch with BLOSUM62, gap_open=-11, gap_extend=-1.
# Same algorithm and parameters as the previous Bio.pairwise2.align.globalds() call.
_ALIGNER = PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.substitution_matrix = SUBSTITUTION_MATRIX
_ALIGNER.open_gap_score = -11
_ALIGNER.extend_gap_score = -1


def _build_blosum_lut(matrix) -> np.ndarray:
    """ASCII-indexed (256, 256) lookup table of BLOSUM62 scores for vectorized scoring."""
    lut = np.zeros((256, 256), dtype=np.int8)
    for a in matrix.alphabet:
        for b in matrix.alphabet:
            lut[ord(a), ord(b)] = matrix[(a, b)]
    return lut


_BLOSUM_LUT = _build_blosum_lut(SUBSTITUTION_MATRIX)
_GAP_BYTE = ord("-")


# --- Parallel execution across CPU cores -------------------------------------
# Bio.Align.PairwiseAligner is GIL-bound (C-impl that does not release the GIL),
# so a ThreadPoolExecutor serializes and runs effectively single-threaded. We use
# *processes* instead to actually saturate the node's cores. The (read-only) train
# set is shared into each worker once via an initializer rather than being pickled
# per task. The aligner / BLOSUM LUT module globals are recreated per worker on
# spawn (or inherited on fork) — either way they exist at module import.
_WORKER_TRAIN: Optional[List[str]] = None
_WORKER_SELF: bool = False

# Pin to a fork context: workers inherit the parent's sys.path + the aligner/LUT
# module globals (no re-import, no per-task repickling of the train set). This is
# the Linux default today but Python plans to change it, so request it explicitly.
# Fall back to the platform default where fork is unavailable (e.g. macOS spawn /
# Windows) — the logic is identical, only the bootstrap differs.
import multiprocessing as _multiprocessing

try:
    _MP_CONTEXT = _multiprocessing.get_context("fork")
except ValueError:  # pragma: no cover - non-fork platforms
    _MP_CONTEXT = None


def _available_cpus() -> int:
    """CPUs actually allocated to this process (respects SLURM/PBS cpusets)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:  # not available on all platforms (e.g. macOS)
        return max(1, os.cpu_count() or 1)


def _init_identity_worker(train: List[str], self_comparison: bool) -> None:
    global _WORKER_TRAIN, _WORKER_SELF
    _WORKER_TRAIN = train
    _WORKER_SELF = self_comparison


def _worker_max(item: Tuple[int, str]):
    """Per-query max identity/similarity against the worker's shared train set."""
    i, generated_seq = item
    local_max_identity = 0.0
    local_max_similarity = 0.0
    local_max_identity_idx = 0
    local_max_similarity_idx = 0
    for j, train_seq in enumerate(_WORKER_TRAIN):
        if _WORKER_SELF and i == j:
            continue
        sequence_identity, sequence_similarity = _pair_metrics(generated_seq, train_seq)
        if sequence_identity > local_max_identity:
            local_max_identity = sequence_identity
            local_max_identity_idx = j
        if sequence_similarity > local_max_similarity:
            local_max_similarity = sequence_similarity
            local_max_similarity_idx = j
    return (
        i,
        local_max_identity,
        local_max_similarity,
        local_max_identity_idx,
        local_max_similarity_idx,
    )


def _worker_topk(item: Tuple[int, str, int]):
    """Per-query top-k reference hits (identity desc) against the shared train set."""
    i, generated_seq, top_k = item
    scored = []
    for j, train_seq in enumerate(_WORKER_TRAIN):
        if _WORKER_SELF and i == j:
            continue
        sequence_identity, _ = _pair_metrics(generated_seq, train_seq)
        scored.append((sequence_identity, j))
    # Sort by identity desc; tie-break on train index asc for determinism.
    scored.sort(key=lambda t: (-t[0], t[1]))
    # Emit identity as PERCENT (native metric for this tool's top-k contract).
    return i, [(j, identity * 100.0) for identity, j in scored[:top_k]]


def evaluate_max_sequence_identity(
    train,
    generated,
    generated_identifiers,
    train_identifiers,
    *,
    save_path: Optional[str] = None,
    return_second_max: bool = False,
    top_k: Optional[int] = None,
    topk_save_path: Optional[str] = None,
):
    """Compute top identity/similarity hits for each generated sequence.

    Args:
        train: Reference amino-acid sequences.
        generated: Sequences evaluated against the reference set.
        generated_identifiers: IDs for generated sequences.
        train_identifiers: IDs for train sequences.
        save_path: Optional path to write CSV output.
        return_second_max: Skip self-comparison for same-index pairs.
        top_k: When >= 1, also write a tidy top-k CSV (columns
            query_id,rank,neighbour_id,score) of the k highest-identity reference
            hits per query. ``score`` is identity PERCENT in [0, 100] (LARGER =
            closer), i.e. 100 * the fraction stored in the default CSV's
            ``sequence_identity`` column.
        topk_save_path: Path for the top-k CSV (required when ``top_k`` is set).
    """
    (
        max_sequence_identity,
        max_sequence_similarity,
        max_sequence_identity_index,
        max_sequence_similarity_index,
    ) = get_max_sequence_identity(
        train,
        generated,
        self_comparison=return_second_max,
    )

    max_sequence_identity_hit = [train_identifiers[idx] for idx in max_sequence_identity_index]
    max_sequence_similarity_hit = [train_identifiers[idx] for idx in max_sequence_similarity_index]

    if save_path is not None:
        df = pd.DataFrame(
            {
                "ID": generated_identifiers,
                "sequence_identity": max_sequence_identity,
                "sequence_identity_hit": max_sequence_identity_hit,
                "sequence_similarity": max_sequence_similarity,
                "sequence_similarity_hit": max_sequence_similarity_hit,
            }
        )
        df.to_csv(save_path, index=False)

    if top_k is not None and top_k >= 1 and topk_save_path is not None:
        topk = get_topk_sequence_identity(
            train,
            generated,
            top_k,
            self_comparison=return_second_max,
        )
        _write_topk(
            topk,
            generated_identifiers,
            train_identifiers,
            topk_save_path,
        )

    return (
        max_sequence_identity,
        max_sequence_similarity,
        max_sequence_identity_hit,
        max_sequence_similarity_hit,
    )


def get_topk_sequence_identity(
    train: Sequence[str],
    generated: Sequence[str],
    top_k: int,
    *,
    self_comparison: bool = False,
):
    """Per generated sequence, the top-k reference hits by identity (descending).

    Returns a list (per query) of lists of (train_index, identity) tuples,
    highest identity first, at most ``top_k`` long. In self mode the query's own
    index is excluded from its neighbour list.
    """
    train = [str(seq) for seq in train]
    generated = [str(seq) for seq in generated]

    results = [[] for _ in generated]
    if not generated:
        return results

    max_workers = max(1, min(len(generated), _available_cpus()))
    items = [(i, seq, top_k) for i, seq in enumerate(generated)]
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=_MP_CONTEXT,
        initializer=_init_identity_worker,
        initargs=(train, self_comparison),
    ) as executor:
        for i, ranked in executor.map(_worker_topk, items, chunksize=8):
            results[i] = ranked

    return results


def _write_topk(topk, generated_identifiers, train_identifiers, save_path: str) -> None:
    rows = []
    for query_idx, ranked in enumerate(topk):
        query_id = str(generated_identifiers[query_idx]).split()[0] if generated_identifiers[query_idx] else generated_identifiers[query_idx]
        for rank, (train_idx, score) in enumerate(ranked, start=1):
            neighbour_id = str(train_identifiers[train_idx]).split()[0]
            rows.append(
                {
                    "query_id": query_id,
                    "rank": rank,
                    "neighbour_id": neighbour_id,
                    "score": score,
                }
            )
    pd.DataFrame(rows, columns=["query_id", "rank", "neighbour_id", "score"]).to_csv(
        save_path, index=False
    )


def _pair_metrics(seq1: str, seq2: str) -> Tuple[float, float]:
    """Return (identity, similarity) for two sequences via global BLOSUM62 alignment.

    Uses Needleman-Wunsch with BLOSUM62, gap_open=-11, gap_extend=-1 — same algorithm
    and parameters as the previous Bio.pairwise2-based implementation. The current
    Bio.Align.PairwiseAligner C-impl is ~18x faster per pair.

    Identity   = (positions where aligned chars match exactly, excluding gaps) / alignment length.
    Similarity = (positions with positive BLOSUM62 score, excluding gaps)      / alignment length.

    Note: when a pair admits multiple equally-optimal alignments, PairwiseAligner may
    pick a different traceback than pairwise2 did, leading to per-pair differences in
    identity/similarity of up to ~2 percentage points (alignment scores are identical).
    Aggregate statistics (mean/median over many pairs) are unaffected to ~5 decimals.
    """
    alignment = _ALIGNER.align(seq1, seq2)[0]
    aligned_seq_1 = str(alignment[0])
    aligned_seq_2 = str(alignment[1])
    denominator = len(aligned_seq_1)

    arr_1 = np.frombuffer(aligned_seq_1.encode("ascii"), dtype=np.uint8)
    arr_2 = np.frombuffer(aligned_seq_2.encode("ascii"), dtype=np.uint8)
    no_gap = (arr_1 != _GAP_BYTE) & (arr_2 != _GAP_BYTE)
    matches = int(((arr_1 == arr_2) & no_gap).sum())
    positives = int(((_BLOSUM_LUT[arr_1, arr_2] > 0) & no_gap).sum())

    return matches / denominator, positives / denominator


def get_max_sequence_identity(
    train,
    generated=None,
    *,
    self_comparison: bool = False,
):
    if generated is None:
        generated = train
        self_comparison = True
    return get_max_sequence_identity_two_sets(
        train,
        generated,
        self_comparison=self_comparison,
    )


def get_max_sequence_identity_two_sets(
    train: Sequence[str],
    generated: Sequence[str],
    *,
    self_comparison: bool = False,
):
    train = [str(seq) for seq in train]
    generated = [str(seq) for seq in generated]

    max_sequence_identity = [0.0 for _ in generated]
    max_sequence_similarity = [0.0 for _ in generated]
    max_sequence_identity_index = [0 for _ in generated]
    max_sequence_similarity_index = [0 for _ in generated]

    if not generated:
        return (
            max_sequence_identity,
            max_sequence_similarity,
            max_sequence_identity_index,
            max_sequence_similarity_index,
        )

    max_workers = max(1, min(len(generated), _available_cpus()))
    items = list(enumerate(generated))
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=_MP_CONTEXT,
        initializer=_init_identity_worker,
        initargs=(train, self_comparison),
    ) as executor:
        for (
            i,
            local_max_identity,
            local_max_similarity,
            local_max_identity_idx,
            local_max_similarity_idx,
        ) in executor.map(_worker_max, items, chunksize=8):
            max_sequence_identity[i] = local_max_identity
            max_sequence_similarity[i] = local_max_similarity
            max_sequence_identity_index[i] = local_max_identity_idx
            max_sequence_similarity_index[i] = local_max_similarity_idx

    return (
        max_sequence_identity,
        max_sequence_similarity,
        max_sequence_identity_index,
        max_sequence_similarity_index,
    )


def _get_save_path(data_path: str, *, save_suffix: Optional[str] = None) -> str:
    extension = data_path.split(".")[-1]
    base_path = data_path[: -len(extension) - 1]
    suffix = "" if save_suffix is None else f"_{save_suffix}"
    return f"{base_path}_max_sequence_identity{suffix}.csv"


def _get_topk_save_path(data_path: str) -> str:
    extension = data_path.split(".")[-1]
    base_path = data_path[: -len(extension) - 1]
    return f"{base_path}_max_sequence_identity_topk.csv"


def max_sequence_identity(
    fasta_path: str,
    *,
    train_path: Optional[str] = None,
    top_k: Optional[int] = None,
):
    if train_path is None:
        main_train_sequences(fasta_path, top_k=top_k)
    else:
        main_generated_sequences(fasta_path=fasta_path, train_path=train_path, top_k=top_k)


def main_generated_sequences(*, fasta_path: str, train_path: str, top_k: Optional[int] = None):
    generated_identifiers, generated = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    train_identifiers, train = separate_identifiers(
        load_fasta_sequences(train_path, load_identifiers=True)
    )
    save_path = _get_save_path(fasta_path)
    evaluate_max_sequence_identity(
        train,
        generated,
        generated_identifiers,
        train_identifiers,
        save_path=save_path,
        top_k=top_k,
        topk_save_path=_get_topk_save_path(fasta_path),
    )


def main_train_sequences(train_path: str, *, top_k: Optional[int] = None):
    train_identifiers, train = separate_identifiers(
        load_fasta_sequences(train_path, load_identifiers=True)
    )
    _main_train_sequences(train_path, train, train_identifiers, top_k=top_k)


def _main_train_sequences(train_path: str, train, train_identifiers, *, top_k: Optional[int] = None):
    save_path = _get_save_path(train_path, save_suffix="self")
    evaluate_max_sequence_identity(
        train,
        train,
        train_identifiers,
        train_identifiers,
        save_path=save_path,
        return_second_max=True,
        top_k=top_k,
        topk_save_path=_get_topk_save_path(train_path),
    )
