"""Unit tests for prepare_order.py (design -> codon-optimized, overhang-flanked DNA).

Run: python test_prepare_order.py

Requires dnachisel (tps_eval env). Uses tiny synthetic FASTA/CSV inputs and hardcoded
DNA. Exercises the pure helpers (revcomp, homopolymer run, GC window, Type IIS overlap
detection), input loading (FASTA + CSV), construct validation, and the full
prepare_one / prepare_order assembly with a translation round-trip.
"""
from __future__ import annotations

import os
import sys
import tempfile

import pandas as pd
from Bio.Seq import Seq


from tps_eval.order_preparation.overhangs import get_overhangs
from tps_eval.order_preparation.prepare_order import (
    _gc_window_extremes,
    _max_homopolymer_run,
    _revcomp,
    _typeiis_violations,
    load_designs,
    normalize_termini,
    prepare_one,
    prepare_order,
    validate_construct,
)


def test_revcomp():
    assert _revcomp("ATGC") == "GCAT"
    assert _revcomp("AAAA") == "TTTT"
    assert _revcomp("atgc") == "gcat"     # case preserved
    print("ok revcomp")


def test_max_homopolymer_run():
    assert _max_homopolymer_run("") == 0
    assert _max_homopolymer_run("ACGT") == 1
    assert _max_homopolymer_run("AATTTG") == 3
    assert _max_homopolymer_run("GGGGGG") == 6
    print("ok max_homopolymer_run")


def test_gc_window_extremes():
    # Shorter than window -> overall GC for both bounds.
    lo, hi = _gc_window_extremes("GCGC", 10)
    assert lo == hi == 1.0
    # 4-bp windows over ATATGCGC: min 0.0 (ATAT), max 1.0 (GCGC).
    lo, hi = _gc_window_extremes("ATATGCGC", 4)
    assert abs(lo - 0.0) < 1e-9 and abs(hi - 1.0) < 1e-9
    print("ok gc_window_extremes")


def test_typeiis_violations():
    prefix, suffix = get_overhangs("Type 3")
    # A clean CDS (no internal site) -> only the deliberate flank sites -> no violation.
    clean_cds = "ATG" + "AAA" * 10 + "TAA"
    full = prefix + clean_cds + suffix
    assert _typeiis_violations(full, prefix, suffix, ("BsaI", "BsmBI")) == []
    # Inject a BsaI site (GGTCTC) inside the CDS -> a violation must be reported.
    bad_cds = "ATG" + "GGTCTC" + "AAA" * 6 + "TAA"
    full_bad = prefix + bad_cds + suffix
    viol = _typeiis_violations(full_bad, prefix, suffix, ("BsaI", "BsmBI"))
    assert viol and any("BsaI" in v for v in viol), viol
    print("ok typeiis_violations")


def test_load_designs_fasta_and_csv():
    tmp = tempfile.mkdtemp(prefix="order_load_")
    fa = os.path.join(tmp, "d.fasta")
    with open(fa, "w") as fh:
        fh.write(">d1\nMAAK\n>d2\nMKLV\n")
    assert load_designs(fa) == [("d1", "MAAK"), ("d2", "MKLV")]

    csv = os.path.join(tmp, "d.csv")
    pd.DataFrame({"id": ["x1", "x2"], "sequence": ["MAAK", "MKLV"]}).to_csv(csv, index=False)
    assert load_designs(csv) == [("x1", "MAAK"), ("x2", "MKLV")]

    # CSV without an id column -> row index used as id.
    csv2 = os.path.join(tmp, "d2.csv")
    pd.DataFrame({"protein": ["MAAK", "MKLV"]}).to_csv(csv2, index=False)
    assert load_designs(csv2) == [("0", "MAAK"), ("1", "MKLV")]
    print("ok load_designs_fasta_and_csv")


def test_normalize_termini():
    # Already a complete ORF at the AA level: begins with Met, ends with an explicit stop.
    assert normalize_termini("MAAK*") == ("MAAK", False, False)
    # Met start, no terminal stop -> a stop is being added (the normal design case).
    assert normalize_termini("MAAK") == ("MAAK", False, True)
    # No Met start, explicit stop -> a start Met is prepended, stop already present.
    assert normalize_termini("AAK*") == ("MAAK", True, False)
    # Neither terminus present -> both added.
    assert normalize_termini("AAK") == ("MAAK", True, True)
    # Whitespace / lowercase are normalized; multiple trailing stops collapse.
    assert normalize_termini("  maak**  ") == ("MAAK", False, False)
    # Empty (or stop-only) input is rejected.
    for bad in ("", "   ", "*", "**"):
        try:
            normalize_termini(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {bad!r}")
    print("ok normalize_termini")


def test_prepare_one_adds_start_and_stop():
    # A design lacking both an N-terminal Met and a terminal stop: prepare_one must add both,
    # flag them, warn, and still round-trip (CDS translates to 'M' + the design).
    design = "AKGEELFTGVVPILVELDGDVNGHK"
    row = prepare_one(design, organism="yeast", overhang_type="Type 3", seed=0)
    assert row["status"] == "ok", row["warnings"]
    assert row["start_added"] is True and row["stop_added"] is True
    assert "start codon" in row["warnings"] and "stop codon" in row["warnings"]
    cds = row["cds"]
    assert cds.upper().startswith("ATG") and cds[-3:].upper() in ("TAA", "TAG", "TGA")
    assert str(Seq(cds).translate(to_stop=False)).rstrip("*") == "M" + design
    print("ok prepare_one_adds_start_and_stop")


def test_prepare_order_reports_added_counts():
    # Mixed batch: one complete-ORF design (Met start, explicit stop) and one missing both.
    # The per-row flags must reflect exactly which terminus each design needed.
    tmp = tempfile.mkdtemp(prefix="order_termini_")
    fa = os.path.join(tmp, "designs.fasta")
    with open(fa, "w") as fh:
        fh.write(">complete\nMASKGEELFTGVV*\n>bare\nAKGEELFTGVV\n")
    df = prepare_order(fa, seed=0, save=False)
    by_id = {r["id"]: r for _, r in df.iterrows()}
    assert by_id["complete"]["start_added"] == False and by_id["complete"]["stop_added"] == False
    assert by_id["bare"]["start_added"] == True and by_id["bare"]["stop_added"] == True
    assert int(df["start_added"].sum()) == 1 and int(df["stop_added"].sum()) == 1
    print("ok prepare_order_reports_added_counts")


def test_prepare_one_assembles_and_roundtrips():
    protein = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGK"
    prefix, suffix = get_overhangs("Type 3")
    row = prepare_one(protein, organism="yeast", overhang_type="Type 3", seed=0)
    assert row["status"] == "ok", row["warnings"]
    full = row["ordered_sequence"]
    assert full.startswith(prefix) and full.endswith(suffix)
    cds = row["cds"]
    assert cds.upper().startswith("ATG") and cds[-3:].upper() in ("TAA", "TAG", "TGA")
    translated = str(Seq(cds).translate(to_stop=False)).rstrip("*")
    assert translated == protein
    # Independent validation must find no hard problems.
    warnings = validate_construct(full, protein, "Type 3")
    hard = [w for w in warnings if "overlaps the CDS" in w or "does not" in w
            or "translated CDS" in w]
    assert not hard, hard
    print("ok prepare_one_assembles_and_roundtrips")


def test_prepare_order_end_to_end():
    tmp = tempfile.mkdtemp(prefix="order_e2e_")
    fa = os.path.join(tmp, "designs.fasta")
    with open(fa, "w") as fh:
        fh.write(">a\nMASKGEELFTGVV\n>b\nMKLVINGKPQDE\n")
    df = prepare_order(fa, seed=0, save=False)
    assert len(df) == 2
    assert (df["status"] == "ok").all(), df[["id", "status", "warnings"]].to_dict()
    for _, r in df.iterrows():
        cds = r["cds"]
        assert str(Seq(cds).translate(to_stop=False)).rstrip("*") == r["protein"]
    print("ok prepare_order_end_to_end")


def test_prepare_order_empty_raises():
    tmp = tempfile.mkdtemp(prefix="order_empty_")
    fa = os.path.join(tmp, "empty.fasta")
    open(fa, "w").close()
    try:
        prepare_order(fa, save=False)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for empty input")
    print("ok prepare_order_empty_raises")


def test_typeiis_reverse_strand_bsmbi():
    # A BsmBI site on the REVERSE strand inside the CDS must be flagged. revcomp(CGTCTC) =
    # GAGACG; note the Type-3 suffix legitimately carries GAGACC/GAGACG (the deliberate
    # reverse-strand flank sites), which must NOT be flagged.
    prefix, suffix = get_overhangs("Type 3")
    clean = prefix + "ATG" + "AAA" * 8 + "TAA" + suffix
    assert _typeiis_violations(clean, prefix, suffix, ("BsaI", "BsmBI")) == []
    bad = prefix + "ATG" + "AAA" * 3 + "GAGACG" + "AAA" * 3 + "TAA" + suffix
    viol = _typeiis_violations(bad, prefix, suffix, ("BsaI", "BsmBI"))
    assert viol and any("BsmBI" in v for v in viol), viol
    print("ok typeiis_reverse_strand_bsmbi")


def test_typeiis_junction_spanning():
    # A site straddling the flank<->CDS junction must be flagged (it would still cut the part).
    # Synthetic flanks make the boundary exact: prefix ends 'GG', CDS starts 'TCTC' ->
    # 'GGTCTC' (BsaI) spans the prefix/CDS junction (cds_start = 6, site at nt 4..10).
    viol = _typeiis_violations("AAAAGG" + "TCTCAAATAA" + "TTTTTT", "AAAAGG", "TTTTTT",
                               ("BsaI", "BsmBI"))
    assert viol and any("BsaI" in v for v in viol), viol
    # The same site sitting WHOLLY inside a flank (ending exactly at cds_start) is allowed.
    assert _typeiis_violations("GGTCTC" + "ATGAAATAA" + "TTTTTT", "GGTCTC", "TTTTTT",
                               ("BsaI", "BsmBI")) == []
    print("ok typeiis_junction_spanning")


def test_validate_construct_flags_internal_site():
    # validate_construct: a clean Type-3 construct -> [] (GC/homopolymer checks disabled so we
    # isolate the Type IIS guarantee); one with an internal BsaI site -> an 'overlaps the CDS'
    # warning. CDSs translate to their stated protein so no translation warning is raised.
    prefix, suffix = get_overhangs("Type 3")
    clean = prefix + "ATGAAGTAA" + suffix                      # translates to 'MK'
    assert validate_construct(clean, "MK", "Type 3",
                              max_homopolymer=0, gc_min=None, gc_max=None) == []
    bad = prefix + "ATGGGTCTCGGGTAA" + suffix                  # 'MGLG', internal BsaI (GGTCTC)
    warnings = validate_construct(bad, "MGLG", "Type 3",
                                  max_homopolymer=0, gc_min=None, gc_max=None)
    assert any("overlaps the CDS" in w for w in warnings), warnings
    print("ok validate_construct_flags_internal_site")


def test_prepare_one_clears_typeiis_site():
    # The codon optimizer (+ its re-optimization loop) must yield a site-free construct even for
    # a site-prone protein (Arg/Gly/Ser/Leu/Pro-rich — codons whose obvious choices include the
    # BsaI/BsmBI motifs). Confirms the guarantee end-to-end, not just the detector.
    prefix, suffix = get_overhangs("Type 3")
    protein = "MARSLRGWLPRSAAGKLERSGPRLGSRL"
    row = prepare_one(protein, organism="yeast", overhang_type="Type 3", seed=0)
    assert row["status"] == "ok", row["warnings"]
    assert _typeiis_violations(row["ordered_sequence"], prefix, suffix,
                               ("BsaI", "BsmBI")) == []
    assert str(Seq(row["cds"]).translate(to_stop=False)).rstrip("*") == protein
    print("ok prepare_one_clears_typeiis_site")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
