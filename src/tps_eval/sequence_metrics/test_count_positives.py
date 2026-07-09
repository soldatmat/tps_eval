from __future__ import annotations

"""Self-contained tests for count_positives.py.

Run from this directory (so the flat-module import resolves like the runner does):
    cd src/sequence_metrics && python test_count_positives.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_count_positives.py -q

count_positives counts aligned residue pairs with a POSITIVE substitution score,
ignoring gap columns. These tests lock in: exact match counting on a known
BLOSUM62 substitution matrix, that gap columns ('-') are skipped on either side,
that a strictly-zero/negative score is NOT counted, and empty/all-gap inputs give
0. Uses only Biopython's bundled BLOSUM62 (no network, no models).
"""

from Bio.Align import substitution_matrices

from tps_eval.sequence_metrics.count_positives import count_positives

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def test_identical_sequence_all_positive():
    # Every diagonal BLOSUM62 score is positive, so all non-gap columns count.
    seq = "ACDEFGHIK"
    assert count_positives(seq, seq, BLOSUM62) == len(seq)


def test_gaps_are_skipped_both_sides():
    # Gap on either strand => that column is ignored entirely.
    a = "A-CD"
    b = "AE-D"
    # Columns: (A,A)=+, (-,E) skip, (C,-) skip, (D,D)=+  => 2 positives.
    assert count_positives(a, b, BLOSUM62) == 2


def test_negative_and_zero_not_counted():
    # A vs W is negative in BLOSUM62 (-3); must not be counted.
    assert BLOSUM62[("A", "W")] < 0
    assert count_positives("A", "W", BLOSUM62) == 0
    # A vs S is +1 (positive) -> counted; construct a mixed pair.
    assert BLOSUM62[("A", "S")] > 0
    assert count_positives("AA", "WS", BLOSUM62) == 1


def test_empty_and_all_gap_zero():
    assert count_positives("", "", BLOSUM62) == 0
    assert count_positives("---", "AAA", BLOSUM62) == 0


def test_zip_truncates_to_shorter():
    # zip() stops at the shorter strand; extra trailing residues are ignored.
    # (aligned strings are normally equal-length; this documents the contract.)
    assert count_positives("AAAA", "AA", BLOSUM62) == 2


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
