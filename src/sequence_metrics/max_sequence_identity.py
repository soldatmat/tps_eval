from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence, Tuple

import os
import pandas as pd
from Bio import pairwise2
from Bio.Align import substitution_matrices

import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers
from count_positives import count_positives


SUBSTITUTION_MATRIX = substitution_matrices.load("BLOSUM62")


def evaluate_max_sequence_identity(
    train,
    generated,
    generated_identifiers,
    train_identifiers,
    *,
    save_path: Optional[str] = None,
    return_second_max: bool = False,
):
    """Compute top identity/similarity hits for each generated sequence.

    Args:
        train: Reference amino-acid sequences.
        generated: Sequences evaluated against the reference set.
        generated_identifiers: IDs for generated sequences.
        train_identifiers: IDs for train sequences.
        save_path: Optional path to write CSV output.
        return_second_max: Skip self-comparison for same-index pairs.
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

    return (
        max_sequence_identity,
        max_sequence_similarity,
        max_sequence_identity_hit,
        max_sequence_similarity_hit,
    )


def _pair_metrics(seq1: str, seq2: str) -> Tuple[float, float]:
    alignment = pairwise2.align.globalds(
        seq1,
        seq2,
        SUBSTITUTION_MATRIX,
        -11,
        -1,
        one_alignment_only=True,
    )[0]

    aligned_seq_1 = alignment.seqA
    aligned_seq_2 = alignment.seqB
    denominator = len(aligned_seq_1)
    matches = sum(
        1
        for a, b in zip(aligned_seq_1, aligned_seq_2)
        if a == b and a != "-" and b != "-"
    )
    positives = count_positives(aligned_seq_1, aligned_seq_2, SUBSTITUTION_MATRIX)
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

    def process_generated(i: int, generated_seq: str):
        local_max_identity = 0.0
        local_max_similarity = 0.0
        local_max_identity_idx = 0
        local_max_similarity_idx = 0

        for j, train_seq in enumerate(train):
            if self_comparison and i == j:
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

    max_workers = max(1, min(len(generated), os.cpu_count() or 1))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for (
            i,
            local_max_identity,
            local_max_similarity,
            local_max_identity_idx,
            local_max_similarity_idx,
        ) in executor.map(lambda item: process_generated(*item), enumerate(generated)):
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


def max_sequence_identity(
    fasta_path: str,
    *,
    train_path: Optional[str] = None,
):
    if train_path is None:
        main_train_sequences(fasta_path)
    else:
        main_generated_sequences(fasta_path=fasta_path, train_path=train_path)


def main_generated_sequences(*, fasta_path: str, train_path: str):
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
    )


def main_train_sequences(train_path: str):
    train_identifiers, train = separate_identifiers(
        load_fasta_sequences(train_path, load_identifiers=True)
    )
    _main_train_sequences(train_path, train, train_identifiers)


def _main_train_sequences(train_path: str, train, train_identifiers):
    save_path = _get_save_path(train_path, save_suffix="self")
    evaluate_max_sequence_identity(
        train,
        train,
        train_identifiers,
        train_identifiers,
        save_path=save_path,
        return_second_max=True,
    )
