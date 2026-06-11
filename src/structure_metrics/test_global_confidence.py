from __future__ import annotations

"""Self-contained tests for global_confidence.py.

Run from this directory (so the flat-module import resolves like the runner does):
    cd src/structure_metrics && python test_global_confidence.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_global_confidence.py -q

These tests write synthetic <ID>_pae.npz files in the shared schema (the SAME schema
esmfold.py / extract_pae.py write) and check the pTM/ipTM extraction + directory
driver in closed form — no fold model needed (numpy + pandas only).
"""

import os
import tempfile

import numpy as np
import pandas as pd

from global_confidence import (
    BASE_COLUMNS,
    IPTM_COLUMN,
    extract_global_confidence_dir,
    load_global_confidence,
)


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _write_npz(path, *, ptm=None, iptm=None, n=4):
    """A shared-schema <ID>_pae.npz; ptm/iptm omitted when None (older-style npz)."""
    pae = np.full((n, n), 5.0, dtype=np.float32)
    fields = dict(
        pae=pae,
        residue_ids=np.arange(1, n + 1, dtype=np.int32),
        n_residues=np.int64(n),
        source="esmfold",
    )
    if ptm is not None:
        fields["ptm"] = np.float32(ptm)
    if iptm is not None:
        fields["iptm"] = np.float32(iptm)
    np.savez_compressed(path, **fields)


def test_load_missing_file_is_nan():
    vals = load_global_confidence("/no/such/x_pae.npz")
    assert np.isnan(vals["ptm"]) and np.isnan(vals[IPTM_COLUMN])


def test_load_ptm_only():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a_pae.npz")
        _write_npz(p, ptm=0.82)
        vals = load_global_confidence(p)
        _approx(vals["ptm"], 0.82)
        assert np.isnan(vals[IPTM_COLUMN])


def test_load_npz_without_ptm_field_is_nan():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "a_pae.npz")
        _write_npz(p)  # no ptm/iptm fields (pre-pTM npz)
        vals = load_global_confidence(p)
        assert np.isnan(vals["ptm"]) and np.isnan(vals[IPTM_COLUMN])


def test_dir_single_chain_no_iptm_column():
    with tempfile.TemporaryDirectory() as d:
        pae_dir = os.path.join(d, "structs_pae")
        os.makedirs(pae_dir)
        _write_npz(os.path.join(pae_dir, "seq1_pae.npz"), ptm=0.7)
        _write_npz(os.path.join(pae_dir, "seq2_pae.npz"), ptm=0.9)
        df = extract_global_confidence_dir(pae_dir)
        assert list(df.columns) == BASE_COLUMNS  # no iptm column
        assert set(df["ID"]) == {"seq1", "seq2"}
        _approx(float(df.set_index("ID").loc["seq1", "ptm"]), 0.7)
        _approx(float(df.set_index("ID").loc["seq2", "ptm"]), 0.9)
        # default name keyed off pae_dir
        assert os.path.isfile(pae_dir + "_global_confidence.csv")


def test_dir_iptm_column_when_present():
    with tempfile.TemporaryDirectory() as d:
        pae_dir = os.path.join(d, "structs_pae")
        os.makedirs(pae_dir)
        _write_npz(os.path.join(pae_dir, "c1_pae.npz"), ptm=0.8, iptm=0.6)
        _write_npz(os.path.join(pae_dir, "c2_pae.npz"), ptm=0.75)  # iptm NaN
        df = extract_global_confidence_dir(pae_dir)
        assert IPTM_COLUMN in df.columns
        _approx(float(df.set_index("ID").loc["c1", IPTM_COLUMN]), 0.6)
        assert pd.isna(df.set_index("ID").loc["c2", IPTM_COLUMN])


def test_dir_named_off_structs_dir_and_fixed_ids():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        pae_dir = os.path.join(d, "structs_pae")
        os.makedirs(structs)
        os.makedirs(pae_dir)
        _write_npz(os.path.join(pae_dir, "have_pae.npz"), ptm=0.65)
        # fixed id set: one present, one missing -> NaN row, every id gets a row
        df = extract_global_confidence_dir(
            pae_dir, structs_dir=structs, ids=["have", "missing"]
        )
        assert set(df["ID"]) == {"have", "missing"}
        _approx(float(df.set_index("ID").loc["have", "ptm"]), 0.65)
        assert pd.isna(df.set_index("ID").loc["missing", "ptm"])
        assert os.path.isfile(structs + "_global_confidence.csv")  # named off structs_dir


def test_save_path_override():
    with tempfile.TemporaryDirectory() as d:
        pae_dir = os.path.join(d, "p")
        os.makedirs(pae_dir)
        _write_npz(os.path.join(pae_dir, "x_pae.npz"), ptm=0.5)
        out = os.path.join(d, "custom.csv")
        extract_global_confidence_dir(pae_dir, save_path=out)
        assert os.path.isfile(out)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
