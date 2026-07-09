from __future__ import annotations

"""Unit tests for the shared TPS-vs-non-TPS accession classification core.

Run: python test_tps_accessions.py   (no pytest dependency required).

Only synthetic data: tiny accession files written to temp dirs and in-memory
frozensets. Locks in accession canonicalization (isoform/version suffix strip),
membership classification, whitespace tolerance, the FileNotFoundError contract,
and that the committed reference set loads.
"""

import os
import sys
import tempfile


from tps_eval.homology_search.tps_accessions import (  # noqa: E402
    _canonical_accession,
    is_tps,
    load_tps_accessions,
)


def _write_accessions(lines) -> str:
    """Write one accession per line to a UNIQUE temp file (load_tps_accessions is
    lru_cached per path, so a fresh path is needed whenever the content changes)."""
    fd, path = tempfile.mkstemp(prefix="tps_acc_", suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def test_canonical_accession():
    assert _canonical_accession("P12345") == "P12345"
    assert _canonical_accession("P12345-2") == "P12345"      # isoform suffix
    assert _canonical_accession("P12345.1") == "P12345"      # version suffix
    assert _canonical_accession("  P12345  ") == "P12345"    # whitespace
    assert _canonical_accession("P12345-2.1") == "P12345"    # isoform then version
    print("ok _canonical_accession")


def test_load_and_is_tps():
    path = _write_accessions(["P0C2A9", "Q9XHX8", "  O48935  "])
    tps = load_tps_accessions(path)
    assert isinstance(tps, frozenset)
    assert tps == {"P0C2A9", "Q9XHX8", "O48935"}, tps      # whitespace trimmed
    # Membership + canonicalization on the query side.
    assert is_tps("P0C2A9", tps) is True
    assert is_tps("P0C2A9-3", tps) is True                 # isoform of a member
    assert is_tps("P0C2A9.2", tps) is True                 # version of a member
    assert is_tps("Q99999", tps) is False                  # non-member
    # A prenyltransferase-style accession that is deliberately NOT a TPS.
    assert is_tps("P99998", tps) is False
    print("ok load + is_tps")


def test_blank_lines_ignored():
    path = _write_accessions(["P0C2A9", "", "   ", "Q9XHX8"])
    tps = load_tps_accessions(path)
    assert tps == {"P0C2A9", "Q9XHX8"}, tps                # blank lines dropped
    print("ok blank lines ignored")


def test_missing_file_raises():
    missing = os.path.join(tempfile.gettempdir(), "definitely_absent_tps_acc_xyz.txt")
    if os.path.exists(missing):
        os.remove(missing)
    try:
        load_tps_accessions(missing)
    except FileNotFoundError:
        print("ok missing file raises FileNotFoundError")
    else:
        raise AssertionError("expected FileNotFoundError for a missing accession file")


def test_committed_reference_set_loads():
    """The committed reference set (src/homology_search/tps_uniprot_accessions.txt)
    loads to a non-trivial frozenset of canonical accessions."""
    here = os.path.dirname(os.path.abspath(__file__))
    ref = os.path.join(here, "tps_uniprot_accessions.txt")
    if os.path.exists(ref):
        tps = load_tps_accessions(ref)
        assert len(tps) > 1000, len(tps)
        # A known reference-set member (from the file head) classifies as TPS.
        assert is_tps("A0A067SEC9", tps) is True
        print(f"ok committed reference set loads ({len(tps)} accessions)")
    else:
        print("skip committed-reference test (file absent)")


def main():
    test_canonical_accession()
    test_load_and_is_tps()
    test_blank_lines_ignored()
    test_missing_file_raises()
    test_committed_reference_set_loads()
    print("\nAll 5 tests passed.")


if __name__ == "__main__":
    main()
