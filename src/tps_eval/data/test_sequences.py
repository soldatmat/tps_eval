from __future__ import annotations

"""Self-contained tests for data/sequences.py (FASTA parsing helpers).

Run from this directory (imports resolve like the runner does):
    cd src/data && python test_sequences.py
or under pytest:
    cd src/data && python -m pytest test_sequences.py -q

Synthetic in-memory FASTA files only. Covers padding removal, id extraction, the
MSA loader (which keeps '-' gaps), and identifier/sequence separation.
"""

import os
import sys
import tempfile
from pathlib import Path


from tps_eval.data.sequences import (  # noqa: E402
    load_fasta_msa,
    load_fasta_sequences,
    read_sequences,
    separate_identifiers,
)


def _write_fasta(rows) -> str:
    tmp = tempfile.mkdtemp(prefix="sequences_")
    path = os.path.join(tmp, "seqs.fasta")
    with open(path, "w") as fh:
        for rid, seq in rows:
            fh.write(f">{rid}\n{seq}\n")
    return path


def test_load_sequences_strips_padding():
    path = _write_fasta([("id1", "ACD-EF"), ("id2", "GH--IK")])
    seqs = load_fasta_sequences(path)          # remove_padding=True by default
    assert seqs == ["ACDEF", "GHIK"]
    print("ok load_sequences_strips_padding")


def test_load_sequences_with_identifiers():
    path = _write_fasta([("id1 extra desc", "ACDEF"), ("id2", "GHIK")])
    recs = load_fasta_sequences(path, load_identifiers=True)
    # Biopython id = first whitespace token.
    assert recs == [("id1", "ACDEF"), ("id2", "GHIK")]
    print("ok load_sequences_with_identifiers")


def test_load_msa_keeps_gaps():
    path = _write_fasta([("id1", "AC-DEF"), ("id2", "AC-DEG")])
    msa = load_fasta_msa(path)
    assert msa == ["AC-DEF", "AC-DEG"]   # alignment gaps preserved
    print("ok load_msa_keeps_gaps")


def test_separate_identifiers():
    ids, seqs = separate_identifiers([("id1", "ACDEF"), ("id2", "GHIK")])
    assert ids == ["id1", "id2"]
    assert seqs == ["ACDEF", "GHIK"]
    print("ok separate_identifiers")


def test_read_sequences_direct():
    from types import SimpleNamespace
    recs = [SimpleNamespace(id="a", seq="AC-D"), SimpleNamespace(id="b", seq="EF")]
    assert read_sequences(recs) == ["ACD", "EF"]
    assert read_sequences(recs, remove_padding=False) == ["AC-D", "EF"]
    assert read_sequences(recs, load_identifiers=True) == [("a", "ACD"), ("b", "EF")]
    print("ok read_sequences_direct")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
