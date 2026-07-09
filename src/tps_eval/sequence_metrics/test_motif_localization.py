from __future__ import annotations

"""Self-contained tests for motif_localization.py — the SHARED source of truth for
the two class-I TPS metal-binding motifs (imported by the sequence and structural
motif tools), so it is tested thoroughly.

Run from this directory (flat-module import resolves like the runner does):
    cd src/sequence_metrics && python test_motif_localization.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_motif_localization.py -q

Locks in: the relaxed DDXXD family regex (D->E substitutions still match), the
NSE/DTE 9-mer regex, 0-based half-open start/end + 1-based start_1 bookkeeping,
FIRST-match (left-to-right) determinism, None on absent/empty, and the
coordinating-residue offset mapping onto the acidic/(S/T)/E residues. All inputs
are synthetic strings.
"""

from tps_eval.sequence_metrics.motif_localization import (
    DDXXD_COORDINATING_OFFSETS,
    DDXXD_PATTERN,
    NSE_DTE_COORDINATING_OFFSETS,
    NSE_DTE_PATTERN,
    MotifMatch,
    coordinating_indices,
    locate_ddxxd,
    locate_nse_dte,
)


def test_ddxxd_strict_match_positions():
    m = locate_ddxxd("DDLLD")
    assert m is not None
    assert m.matched == "DDLLD"
    assert (m.start, m.end, m.start_1) == (0, 5, 1)
    # start/end are a half-open 0-based span: seq[start:end] == matched.
    assert "DDLLD"[m.start:m.end] == m.matched


def test_ddxxd_relaxed_family_de_substitutions():
    # DExxD / EDxxD / EExxE conservative substitutions must still localize.
    for seq in ("DEAAD", "EDAAD", "EEAAE", "DDAAE"):
        m = locate_ddxxd(seq)
        assert m is not None and m.start == 0, seq
        assert m.matched == seq


def test_ddxxd_offset_into_sequence():
    m = locate_ddxxd("GGGDDLLD")
    assert m is not None
    assert (m.start, m.end, m.start_1) == (3, 8, 4)


def test_ddxxd_first_match_wins():
    # Two DDXXD-family motifs; the FIRST (left-most) is the deterministic pick.
    m = locate_ddxxd("DDAAD" + "KKKK" + "EEAAE")
    assert m is not None
    assert m.start == 0 and m.matched == "DDAAD"


def test_ddxxd_absent_and_empty():
    assert locate_ddxxd("GGGGGGGG") is None
    assert locate_ddxxd("") is None


def test_nse_dte_match_positions():
    m = locate_nse_dte("NDLASACDE")
    assert m is not None
    assert m.matched == "NDLASACDE"
    assert (m.start, m.end, m.start_1) == (0, 9, 1)
    # variant with D-start and (I|V) at pos 3, T at the (S/T) slot
    m2 = locate_nse_dte("DDIATACDE")
    assert m2 is not None and m2.matched == "DDIATACDE"


def test_nse_dte_absent_and_empty():
    assert locate_nse_dte("GGGGGGGGGG") is None
    assert locate_nse_dte("") is None


def test_coordinating_offset_constants():
    # These offsets are the contract other tools rely on; pin them.
    assert DDXXD_COORDINATING_OFFSETS == (0, 1, 4)
    assert NSE_DTE_COORDINATING_OFFSETS == (0, 1, 4, 8)


def test_coordinating_indices_land_on_expected_residues():
    seq = "GGGDDLLDXXXXNDLASACDE"
    d = locate_ddxxd(seq)
    n = locate_nse_dte(seq)
    d_idx = coordinating_indices(d, DDXXD_COORDINATING_OFFSETS)
    n_idx = coordinating_indices(n, NSE_DTE_COORDINATING_OFFSETS)
    # DDXXD coordinating residues are the three acidic ones (D/E).
    assert [seq[i] for i in d_idx] == ["D", "D", "D"]
    # NSE/DTE coordinating residues are (N/D), D, (S/T), E.
    assert [seq[i] for i in n_idx] == ["N", "D", "S", "E"]
    # Indices are absolute into the full sequence, in the motif's frame.
    assert d_idx == [d.start + o for o in DDXXD_COORDINATING_OFFSETS]


def test_coordinating_indices_drops_out_of_span_offsets():
    m = MotifMatch(matched="DDLLD", start=10, end=15, start_1=11)
    # Offset 99 is outside the 5-mer span -> dropped; valid ones survive.
    assert coordinating_indices(m, (0, 4, 99)) == [10, 14]


def test_patterns_are_the_documented_regexes():
    assert DDXXD_PATTERN.pattern == r"[DE][DE]..[DE]"
    assert NSE_DTE_PATTERN.pattern == r"(N|D)D(L|I|V).(S|T)...E"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
