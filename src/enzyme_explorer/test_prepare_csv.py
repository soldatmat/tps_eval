"""Unit tests for prepare_csv.py (FASTA -> ID,sequence CSV).

Run: python test_prepare_csv.py   (no pytest dependency required).
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from prepare_csv import fasta_to_csv  # noqa: E402


def test_fasta_to_csv_id_and_sequence():
    tmp = tempfile.mkdtemp(prefix="prepcsv_")
    fa = os.path.join(tmp, "seqs.fasta")
    csv = os.path.join(tmp, "seqs.csv")
    with open(fa, "w") as fh:
        fh.write(">design_1 some description\nMAAK\n>design_2\nMKLV\n")
    fasta_to_csv(fa, csv)
    df = pd.read_csv(csv)
    assert list(df.columns) == ["ID", "sequence"]
    # Biopython uses the first whitespace token as the record id.
    assert list(df["ID"]) == ["design_1", "design_2"]
    assert list(df["sequence"]) == ["MAAK", "MKLV"]
    print("ok fasta_to_csv_id_and_sequence")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
