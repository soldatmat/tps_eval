"""Unit tests for codon_optimization.py (protein -> synthesis-ready CDS).

Run: python test_codon_optimization.py

Requires dnachisel (installed in the tps_eval env). Runs the optimizer FOR REAL on tiny
synthetic peptides and asserts the correctness-critical guarantees: in-frame CDS that
translates back to the input AA, starts with ATG, ends in a stop codon, and carries no
internal BsaI/BsmBI site on either strand.
"""
from __future__ import annotations

import os
import sys

from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from codon_optimization import (  # noqa: E402
    _relaxation_ladder,
    codon_optimize,
    resolve_organism,
)

_STOP = ("TAA", "TAG", "TGA")


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _assert_valid_cds(cds: str, protein: str, add_stop: bool = True):
    cds = cds.upper()
    assert cds.startswith("ATG"), f"CDS must start with ATG, got {cds[:3]}"
    assert len(cds) % 3 == 0, f"CDS length {len(cds)} not a multiple of 3"
    # No internal Golden Gate Type IIS site (either strand).
    for site in ("GGTCTC", "CGTCTC"):
        for pat in (site, _revcomp(site)):
            assert pat not in cds, f"forbidden site {pat} present in CDS"
    if add_stop:
        assert cds[-3:] in _STOP, f"CDS must end in a stop, got {cds[-3:]}"
        translated = str(Seq(cds).translate(to_stop=False)).rstrip("*")
    else:
        translated = str(Seq(cds).translate(to_stop=False))
    assert translated == protein.upper().rstrip("*"), (
        f"round-trip mismatch: {translated!r} != {protein!r}"
    )


def test_resolve_organism():
    assert resolve_organism("yeast") == "s_cerevisiae_4932"
    assert resolve_organism("  Yeast ") == "s_cerevisiae_4932"
    assert resolve_organism("s_cerevisiae") == "s_cerevisiae_4932"
    # Unknown names pass through unchanged (already an identifier/taxid).
    assert resolve_organism("s_cerevisiae_4932") == "s_cerevisiae_4932"
    print("ok resolve_organism")


def test_roundtrip_and_no_enzyme_sites():
    protein = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK"
    cds = codon_optimize(protein, organism="yeast", seed=0)
    _assert_valid_cds(cds, protein)
    print("ok roundtrip_and_no_enzyme_sites")


def test_deterministic_with_seed():
    protein = "MAAKLLDEFGHIKMNPQRSTVWY"
    a = codon_optimize(protein, seed=0)
    b = codon_optimize(protein, seed=0)
    assert a == b, "same seed must give identical CDS"
    _assert_valid_cds(a, protein)
    print("ok deterministic_with_seed")


def test_trailing_star_stripped_and_no_stop_option():
    protein = "MAAK*"
    cds = codon_optimize(protein, seed=0, add_stop=False)
    # add_stop=False -> translation equals the star-stripped protein (no stop stripped).
    _assert_valid_cds(cds, "MAAK", add_stop=False)
    assert len(cds) == len("MAAK") * 3   # no extra stop codon appended
    print("ok trailing_star_stripped_and_no_stop_option")


def test_empty_protein_raises():
    try:
        codon_optimize("   ")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError on empty protein")
    print("ok empty_protein_raises")


def test_relaxation_ladder_starts_strict():
    ladder = _relaxation_ladder(max_homopolymer=6, gc_min=0.30, gc_max=0.65, gc_window=50)
    # First rung is the unrelaxed (note=None) configuration.
    assert ladder[0][0] is None
    assert ladder[0][1:] == (6, 0.30, 0.65, 50)
    # Later rungs carry a human-readable relaxation note and loosen constraints.
    assert any(rung[0] is not None for rung in ladder[1:])
    print("ok relaxation_ladder_starts_strict")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
