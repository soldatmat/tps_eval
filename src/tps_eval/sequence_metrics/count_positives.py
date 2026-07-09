from __future__ import annotations


def count_positives(aligned_seq_a: str, aligned_seq_b: str, substitution_matrix) -> int:
    """Count aligned residue pairs with positive substitution score.

    Gap positions are ignored to mirror the original Julia implementation.
    """
    positive_count = 0
    for a, b in zip(aligned_seq_a, aligned_seq_b):
        if a == "-" or b == "-":
            continue
        if substitution_matrix[(a, b)] > 0:
            positive_count += 1
    return positive_count
