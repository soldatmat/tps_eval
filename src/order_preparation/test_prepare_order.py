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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from overhangs import get_overhangs  # noqa: E402
from prepare_order import (  # noqa: E402
    _gc_window_extremes,
    _max_homopolymer_run,
    _revcomp,
    _typeiis_violations,
    load_designs,
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


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
