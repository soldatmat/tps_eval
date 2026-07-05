"""Unit tests for overhangs.py (Golden Gate / MoClo overhang wrapping).

Run: python test_overhangs.py   (no pytest dependency required).

Locks in the validated Type-3 flanks (cross-checked against a real previous order),
the leading-stop stripping of the right fusion, and the unknown-type guard.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from overhangs import (  # noqa: E402
    DEFAULT_OVERHANG,
    OVERHANGS,
    _strip_leading_stop,
    add_overhangs,
    get_overhangs,
)


def test_type3_exact_flanks():
    # These are the exact flanks cross-checked against all 32 sequences of a real order.
    prefix, suffix = get_overhangs("Type 3")
    assert prefix == "actcgacaacCGTCTCatcGGTCTCaT", prefix
    assert suffix == "ATCCtGAGACCtGAGACGgttgtggtgt", suffix
    # Type 3 is the default.
    assert get_overhangs() == get_overhangs(DEFAULT_OVERHANG) == (prefix, suffix)
    print("ok type3_exact_flanks")


def test_strip_leading_stop():
    # Type-3 right fusion bakes a TAG stop that must be stripped (the CDS supplies the stop).
    assert _strip_leading_stop("TAGATCC") == "ATCC"
    assert _strip_leading_stop("TAAxyz") == "xyz"
    assert _strip_leading_stop("TGAggg") == "ggg"
    # A non-stop-prefixed fusion is untouched.
    assert _strip_leading_stop("AACG") == "AACG"
    assert _strip_leading_stop("TATG") == "TATG"   # TAT is not a stop codon
    print("ok strip_leading_stop")


def test_add_overhangs_wraps_cds():
    cds = "ATGAAATAA"
    prefix, suffix = get_overhangs("Type 3")
    assert add_overhangs(cds, "Type 3") == prefix + cds + suffix
    print("ok add_overhangs_wraps_cds")


def test_unknown_type_raises():
    try:
        get_overhangs("Type 999")
    except KeyError as exc:
        assert "Type 999" in str(exc)
    else:
        raise AssertionError("expected KeyError for unknown overhang type")
    print("ok unknown_type_raises")


def test_registry_keys_match_names():
    # Every registry key must match its OverhangType.name (guards copy/paste drift).
    for key, oh in OVERHANGS.items():
        assert oh.name == key, (key, oh.name)
    print("ok registry_keys_match_names")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
