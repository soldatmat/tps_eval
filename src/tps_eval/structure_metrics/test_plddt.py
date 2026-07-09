from __future__ import annotations

"""Self-contained tests for plddt.py — the CANONICAL structure loader (af_output-vs-
flat auto-detection, ID = filename stem, ``<structs_dir>_plddt.csv`` naming) that the
other structure tools reuse. numpy + pandas + biopython only (no torch/pymol).

Run from this directory so the flat-module imports resolve like the runner does:
    cd src/structure_metrics && python test_plddt.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_plddt.py -q

pLDDT lives in the B-factor field, so the I/O tests write tiny PDB/mmCIF files whose
B-factor holds known per-residue confidence values and assert the summary stats, the
HETATM/altloc/multi-chain handling, the af3-vs-flat layout detection + ID keying, the
sibling-CSV filename, and the NaN-on-broken-structure contract.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from tps_eval.structure_metrics.plddt import (
    COLUMNS,
    CONFIDENT_THRESHOLD,
    _collect_structures,
    _default_save_path,
    extract_plddt_dir,
    residue_plddts,
    summarize,
)


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _atom_line(serial, name, resname, chain, resseq, xyz, *, record="ATOM",
               element=None, bfac=80.0, altloc=" ", occ=1.0):
    """One PDB ATOM/HETATM record with an explicit B-factor (== pLDDT here)."""
    x, y, z = xyz
    if element is None:
        element = name[0]
    atom_field = name[:4] if len(name) >= 4 else " " + name.ljust(3)
    return (
        f"{record:<6}{serial:>5} {atom_field}{altloc:1}{resname:>3} {chain}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}{'':10}{element:>2}\n"
    )


def _write_ca_pdb(path, plddts, *, chain="A", start=1):
    """Write a PDB with one CA (ALA) per residue whose B-factor is the given pLDDT."""
    with open(path, "w") as fh:
        for i, b in enumerate(plddts):
            fh.write(_atom_line(i + 1, "CA", "ALA", chain, start + i,
                                (i * 3.8, 0.0, 0.0), element="C", bfac=b))
        fh.write("END\n")


# --- minimal AF3-style mmCIF (B_iso holds pLDDT) ------------------------------- #
_CIF_HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
"""


def _write_ca_cif(path, plddts, *, chain="A"):
    lines = [_CIF_HEADER]
    for i, b in enumerate(plddts):
        lines.append(
            f"ATOM {i+1} C CA . ALA {chain} {i+1} ? {i*3.8:.3f} 0.0 0.0 1.0 {b:.2f} {i+1} {chain} 1\n"
        )
    with open(path, "w") as fh:
        fh.writelines(lines)


# ------------------------------ summarize ------------------------------------- #
def test_summarize_known_values():
    m = summarize([90.0, 80.0, 70.0, 40.0])
    _approx(m["mean_plddt"], 70.0)
    _approx(m["median_plddt"], 75.0)
    _approx(m["min_plddt"], 40.0)
    # 3 of 4 residues are >= 70 (the confident threshold).
    _approx(m["frac_plddt_confident"], 0.75)
    assert m["n_residues"] == 4


def test_summarize_empty_is_nan():
    m = summarize([])
    for k in ("mean_plddt", "median_plddt", "min_plddt", "frac_plddt_confident"):
        assert np.isnan(m[k]), k
    assert m["n_residues"] == 0


def test_summarize_threshold_boundary():
    # >= threshold counts as confident (inclusive).
    m = summarize([CONFIDENT_THRESHOLD, CONFIDENT_THRESHOLD - 0.01])
    _approx(m["frac_plddt_confident"], 0.5)
    m2 = summarize([50.0, 60.0], confident_threshold=55.0)
    _approx(m2["frac_plddt_confident"], 0.5)


# ------------------------------ residue_plddts -------------------------------- #
def test_residue_plddts_reads_ca_bfactor_and_skips_hetatm():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "design.pdb")
        with open(p, "w") as fh:
            fh.write(_atom_line(1, "CA", "ALA", "A", 1, (0, 0, 0), element="C", bfac=95.0))
            fh.write(_atom_line(2, "CA", "GLY", "A", 2, (4, 0, 0), element="C", bfac=55.0))
            # A HETATM calcium ion whose atom is literally named CA must be ignored
            # (else its B-factor would be miscounted as a residue pLDDT).
            fh.write(_atom_line(999, "CA", "CA", "A", 900, (99, 99, 99),
                                record="HETATM", element="CA", bfac=12.0))
            fh.write("END\n")
        got = residue_plddts(p)
        assert got == [95.0, 55.0], got


def test_altloc_ca_counted_once():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "alt.pdb")
        with open(p, "w") as fh:
            # Two altlocs of the same CA: must count ONCE (higher-occupancy wins).
            fh.write(_atom_line(1, "CA", "ALA", "A", 1, (0, 0, 0), element="C",
                                bfac=88.0, altloc="A", occ=0.6))
            fh.write(_atom_line(2, "CA", "ALA", "A", 1, (0.1, 0, 0), element="C",
                                bfac=50.0, altloc="B", occ=0.4))
            fh.write(_atom_line(3, "CA", "GLY", "A", 2, (4, 0, 0), element="C", bfac=60.0))
            fh.write("END\n")
        got = residue_plddts(p)
        assert got == [88.0, 60.0], got


def test_residue_plddts_multichain():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "multi.pdb")
        with open(p, "w") as fh:
            fh.write(_atom_line(1, "CA", "ALA", "A", 1, (0, 0, 0), element="C", bfac=90.0))
            fh.write(_atom_line(2, "CA", "ALA", "B", 1, (0, 5, 0), element="C", bfac=70.0))
            fh.write(_atom_line(3, "CA", "ALA", "B", 2, (4, 5, 0), element="C", bfac=50.0))
            fh.write("END\n")
        got = residue_plddts(p)
        assert sorted(got) == [50.0, 70.0, 90.0], got


def test_residue_plddts_cif():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.cif")
        _write_ca_cif(p, [88.0, 60.0, 30.0])
        got = residue_plddts(p)
        assert got == [88.0, 60.0, 30.0], got


# ------------------------------ _collect_structures --------------------------- #
def test_collect_flat_pdb_wins_over_cif():
    with tempfile.TemporaryDirectory() as d:
        _write_ca_pdb(os.path.join(d, "a.pdb"), [90.0])
        _write_ca_cif(os.path.join(d, "a.cif"), [10.0])   # same stem, .pdb should win
        _write_ca_cif(os.path.join(d, "b.cif"), [70.0])
        structures, mode = _collect_structures(d)
        assert mode == "flat"
        assert list(structures.keys()) == ["a", "b"]
        assert structures["a"].endswith("a.pdb")
        assert structures["b"].endswith("b.cif")


def test_collect_af3_layout_takes_precedence():
    with tempfile.TemporaryDirectory() as d:
        # AF3 af_output: one subfolder per job, each with <job>/<job>_model.cif.
        for job in ("seq_0", "seq_1"):
            sub = os.path.join(d, job)
            os.makedirs(sub)
            _write_ca_cif(os.path.join(sub, job + "_model.cif"), [80.0, 80.0])
        # A stray flat .pdb must be ignored once the af3 layout is detected.
        _write_ca_pdb(os.path.join(d, "decoy.pdb"), [1.0])
        structures, mode = _collect_structures(d)
        assert mode == "af3", mode
        assert list(structures.keys()) == ["seq_0", "seq_1"]
        assert structures["seq_0"].endswith(os.path.join("seq_0", "seq_0_model.cif"))


def test_collect_missing_dir_is_empty():
    structures, mode = _collect_structures("/no/such/dir/hopefully")
    assert len(structures) == 0
    assert mode == "flat"


def test_default_save_path_sibling_and_trailing_slash():
    assert _default_save_path("/x/y/structs") == "/x/y/structs_plddt.csv"
    # A trailing separator must not change the sibling CSV path.
    assert _default_save_path("/x/y/structs/") == "/x/y/structs_plddt.csv"


# ------------------------------ extract_plddt_dir ----------------------------- #
def test_extract_dir_id_keying_csv_and_nan_on_broken():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        _write_ca_pdb(os.path.join(structs, "good.pdb"), [90.0, 80.0, 40.0])
        # Broken: no parseable ATOM records -> NaN row, must NOT abort the batch.
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("this is not a pdb\n")
        df = extract_plddt_dir(structs)

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        # Default save path is the sibling CSV.
        assert os.path.isfile(structs + "_plddt.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["mean_plddt"]), 70.0)
        _approx(float(good["min_plddt"]), 40.0)
        assert int(good["n_residues"]) == 3

        broken = df.set_index("ID").loc["broken"]
        assert np.isnan(broken["mean_plddt"])
        assert int(broken["n_residues"]) == 0


def test_extract_dir_af3_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "af_output")
        os.makedirs(structs)
        for job, plddts in (("seq_0", [90.0, 90.0]), ("seq_1", [40.0, 60.0])):
            sub = os.path.join(structs, job)
            os.makedirs(sub)
            _write_ca_cif(os.path.join(sub, job + "_model.cif"), plddts)
        out = os.path.join(d, "out.csv")
        df = extract_plddt_dir(structs, save_path=out)
        assert os.path.isfile(out)
        assert list(df["ID"]) == ["seq_0", "seq_1"]
        _approx(float(df.set_index("ID").loc["seq_0"]["mean_plddt"]), 90.0)


def test_extract_dir_empty_raises():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "empty")
        os.makedirs(structs)
        try:
            extract_plddt_dir(structs)
        except ValueError:
            return
        raise AssertionError("expected ValueError on an empty structures dir")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
