from __future__ import annotations

"""Self-contained tests for build_cofold_input.py (the AF3 fan-out input builder).

Run from this directory (so the package-style imports inside the module resolve):
    cd src/alphafold && python test_build_cofold_input.py
or under pytest:
    cd src/alphafold && python -m pytest test_build_cofold_input.py -q

Pure csv/tempfile I/O — no AF3, no GPU, no network. The mg_ee grouping is exercised
WITHOUT the EnzymeExplorer/knn dependency by monkeypatching the module-level
`_ee_substrate_per_design` (the seq-only EE argmax loader) with a hardcoded dict, so
the grouping/fallback/manifest logic is tested in isolation. Locks in: ion sets per
mode, one CSV + manifest for the single-group modes, correct headers/rows, the mg_ee
group-per-substrate fan-out + Mg-only fallback, and manifest row correctness.
"""

import csv
import os
import sys
import tempfile
from pathlib import Path


import tps_eval.alphafold.build_cofold_input as bci
from tps_eval.alphafold.build_cofold_input import build, ions_for, read_fasta
from tps_eval.alphafold.cofold_substrates import SUBSTRATE_SMILES


def _write_fasta(path, recs):
    with open(path, "w") as fh:
        for rid, seq in recs:
            fh.write(f">{rid}\n{seq}\n")


def _read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh))


def _read_manifest(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split("\t")
        assert header == ["csv_path", "has_ligand", "n_designs"], header
        for line in fh:
            p, has_lig, n = line.rstrip("\n").split("\t")
            rows.append((p, int(has_lig), int(n)))
    return rows


def test_ions_for_modes():
    assert ions_for("none") == []
    assert ions_for("mg") == [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG")]
    assert ions_for("mg_ppi") == [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG"), ("PPI", "POP")]
    # A forced-substrate mode carries only the 3 Mg (its diphosphate rides on the ligand).
    assert ions_for("mg_gpp") == [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG")]


def test_read_fasta_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        with open(fa, "w") as fh:
            # multi-line sequence + a header with a description (only first token kept)
            fh.write(">d1 some description\nMKT\nAAR\n>d2\nGGG\n")
        recs = read_fasta(fa)
        assert recs == [("d1", "MKTAAR"), ("d2", "GGG")], recs


def test_none_mode_no_ions_no_ligand():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK"), ("b", "GG")])
        manifest = build(fa, "none", d)
        assert manifest == [(os.path.join(d, "af3_input.csv"), False, 2)]
        rows = _read_csv(manifest[0][0])
        assert rows[0] == ["ID", "sequence"], rows[0]
        assert rows[1] == ["a", "MK"]
        assert rows[2] == ["b", "GG"]
        man = _read_manifest(os.path.join(d, "af3_cofold_manifest.tsv"))
        assert man == [(manifest[0][0], 0, 2)]


def test_mg_mode_ion_columns():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK")])
        manifest = build(fa, "mg", d)
        rows = _read_csv(manifest[0][0])
        assert rows[0] == ["ID", "sequence", "ion1_id", "ion1_ccd",
                           "ion2_id", "ion2_ccd", "ion3_id", "ion3_ccd"], rows[0]
        assert rows[1] == ["a", "MK", "MG1", "MG", "MG2", "MG", "MG3", "MG"]
        assert manifest[0][1] is False  # no ligand column


def test_mg_ppi_adds_pop_ion_no_ligand():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK")])
        manifest = build(fa, "mg_ppi", d)
        header = _read_csv(manifest[0][0])[0]
        assert "ion4_id" in header and "ion4_ccd" in header
        row = _read_csv(manifest[0][0])[1]
        # POP placed as the 4th ion (CCD), not as a ligand.
        assert row == ["a", "MK", "MG1", "MG", "MG2", "MG", "MG3", "MG", "PPI", "POP"]
        assert "lig1_smiles" not in header
        assert manifest[0][1] is False


def test_forced_substrate_mode_adds_ligand():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK"), ("b", "GG")])
        manifest = build(fa, "mg_fpp", d)
        assert manifest == [(os.path.join(d, "af3_input.csv"), True, 2)]
        rows = _read_csv(manifest[0][0])
        assert rows[0][-2:] == ["lig1_id", "lig1_smiles"], rows[0]
        # every design gets the SAME forced substrate SMILES
        assert rows[1][-2:] == ["LIG", SUBSTRATE_SMILES["FPP"]]
        assert rows[2][-2:] == ["LIG", SUBSTRATE_SMILES["FPP"]]
        man = _read_manifest(os.path.join(d, "af3_cofold_manifest.tsv"))
        assert man == [(manifest[0][0], 1, 2)]


def test_unknown_mode_raises():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK")])
        try:
            build(fa, "mg_bogus", d)
        except SystemExit:
            return
        raise AssertionError("expected SystemExit for unknown mode")


def test_mg_ee_requires_csv():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("a", "MK")])
        try:
            build(fa, "mg_ee", d, ee_csv=None)
        except SystemExit:
            return
        raise AssertionError("expected SystemExit when mg_ee has no --enzymeexplorer_csv")


def test_mg_ee_groups_by_substrate_with_mgonly_fallback(monkeypatch=None):
    """mg_ee fans out into one CSV per co-foldable substrate + a Mg-only group for
    designs whose EE argmax is not co-foldable. Manifest must carry one row per group
    with the right has_ligand flag and design count."""
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [("g1", "AAA"), ("g2", "CCC"), ("f1", "DDD"), ("x1", "EEE")])
        # Two GPP designs, one FPP, one non-cofoldable (EDSQ) -> Mg-only.
        fake = {"g1": "GPP", "g2": "GPP", "f1": "FPP", "x1": "EDSQ"}
        orig = bci._ee_substrate_per_design
        bci._ee_substrate_per_design = lambda ee_csv: fake
        try:
            manifest = build(fa, "mg_ee", d, ee_csv=os.path.join(d, "ee.csv"))
        finally:
            bci._ee_substrate_per_design = orig

        by_name = {os.path.basename(p): (has, n) for p, has, n in manifest}
        assert set(by_name) == {"af3_input_gpp.csv", "af3_input_fpp.csv", "af3_input_mgonly.csv"}, by_name
        assert by_name["af3_input_gpp.csv"] == (True, 2)
        assert by_name["af3_input_fpp.csv"] == (True, 1)
        assert by_name["af3_input_mgonly.csv"] == (False, 1)

        # GPP group carries the GPP SMILES for both designs.
        gpp_rows = _read_csv(os.path.join(d, "af3_input_gpp.csv"))
        assert gpp_rows[0][-2:] == ["lig1_id", "lig1_smiles"]
        assert {r[0] for r in gpp_rows[1:]} == {"g1", "g2"}
        for r in gpp_rows[1:]:
            assert r[-1] == SUBSTRATE_SMILES["GPP"]

        # Mg-only group has NO ligand column and holds the non-cofoldable design.
        mg_rows = _read_csv(os.path.join(d, "af3_input_mgonly.csv"))
        assert "lig1_smiles" not in mg_rows[0]
        assert [r[0] for r in mg_rows[1:]] == ["x1"]

        # Manifest on disk mirrors the returned manifest (one row per group).
        man = _read_manifest(os.path.join(d, "af3_cofold_manifest.tsv"))
        assert len(man) == 3
        assert sum(n for _, _, n in man) == 4  # every design accounted for exactly once


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
