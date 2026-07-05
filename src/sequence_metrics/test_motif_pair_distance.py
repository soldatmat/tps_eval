from __future__ import annotations

"""Self-contained tests for motif_pair_distance.py.

Run from this directory (flat-module import resolves like the runner does; the
module inserts src/ onto sys.path for its `from data.sequences` import):
    cd src/sequence_metrics && python test_motif_pair_distance.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_motif_pair_distance.py -q

Locks in the residue-separation math between the DDXXD and NSE/DTE motifs
(`motif_start_distance` = 1-based start separation; `motif_gap` = 0-based
inter-motif gap), the NaN-when-a-motif-is-absent contract (a bad row does NOT
abort the run), ID keying, the `<input>_motif_pair_distance.csv` save path, and
the empty-FASTA robustness (regression guard for `pd.DataFrame(rows)[COLUMNS]`
crashing on zero rows). Synthetic FASTA strings only.
"""

import os
import tempfile

import numpy as np

from motif_pair_distance import COLUMNS, _row, motif_pair_distance


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for ident, seq in records:
            fh.write(f">{ident}\n{seq}\n")


# Canonical layout: DDXXD at 0-based [0,5), NSE/DTE at [10,19).
BOTH = "DDLLD" + "GGGGG" + "NDLASACDE"


def test_row_both_motifs_distance_math():
    r = _row("x", BOTH)
    assert r["ddxxd_motif"] == "DDLLD" and r["ddxxd_start"] == 1
    assert r["nse_dte_motif"] == "NDLASACDE" and r["nse_dte_start"] == 11
    # start separation (1-based): 11 - 1 = 10
    assert r["motif_start_distance"] == 10
    # inter-motif gap (0-based): start_NSE(10) - end_DDXXD(5) = 5
    assert r["motif_gap"] == 5


def test_row_missing_nse_gives_nan_distances():
    r = _row("x", "DDLLD" + "GGGGGGGG")  # DDXXD only
    assert r["ddxxd_motif"] == "DDLLD"
    assert r["nse_dte_motif"] == ""
    assert np.isnan(r["motif_start_distance"])
    assert np.isnan(r["motif_gap"])
    assert np.isnan(r["nse_dte_start"])


def test_row_no_motifs_all_nan():
    r = _row("x", "GGGGGGGGGG")
    assert r["ddxxd_motif"] == "" and r["nse_dte_motif"] == ""
    assert np.isnan(r["motif_start_distance"]) and np.isnan(r["motif_gap"])
    assert np.isnan(r["ddxxd_start"]) and np.isnan(r["nse_dte_start"])


def test_dataframe_keyed_by_id_and_saved():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "seqs.fasta")
        _write_fasta(p, [("good", BOTH), ("bad", "GGGGGGGG")])
        df = motif_pair_distance(p, save=True)
        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "bad"}
        by_id = df.set_index("ID")
        assert by_id.loc["good", "motif_start_distance"] == 10
        # Malformed/absent-motif row is NaN, not a crash.
        assert np.isnan(by_id.loc["bad", "motif_start_distance"])
        assert os.path.isfile(os.path.join(d, "seqs_motif_pair_distance.csv"))


def test_id_uses_first_token():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "seqs.fasta")
        _write_fasta(p, [("seq1 description here", BOTH)])
        df = motif_pair_distance(p, save=False)
        assert list(df["ID"]) == ["seq1"]


def test_empty_fasta_returns_empty_frame_with_columns():
    # Regression: zero-sequence FASTA must yield an empty COLUMNS frame, not raise.
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "empty.fasta")
        open(p, "w").close()
        df = motif_pair_distance(p, save=True)
        assert list(df.columns) == COLUMNS
        assert len(df) == 0
        assert os.path.isfile(os.path.join(d, "empty_motif_pair_distance.csv"))


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
