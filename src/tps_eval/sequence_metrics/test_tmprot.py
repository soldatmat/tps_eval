from __future__ import annotations

"""Self-contained tests for tmprot.py (TmProt Tm reshaping — PURE-PYTHON parts).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/sequence_metrics && python test_tmprot.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_tmprot.py -q

The `tmprot` CLI (ESM2 + torch) is NEVER invoked. We monkeypatch subprocess.run to
(a) capture the command that WOULD be launched and (b) drop a synthetic
<out_dir>/<stem>.csv with TmProt's native columns, then assert:
  * run_tmprot_cli builds the right argv (-i/-o/-d ,) and returns the output path,
  * a nonzero exit and a missing output CSV both raise RuntimeError,
  * score_fasta keys by ID, keeps COLUMNS order, sorts by ID, names the CSV
    <fasta_stem>_tmprot.csv, drops Rank/Thermostable, and emits NaN for FASTA
    IDs that TmProt skipped (absent from its output).
"""

import os
import tempfile

import numpy as np
import pandas as pd

import tps_eval.sequence_metrics.tmprot as tmprot
from tps_eval.sequence_metrics.tmprot import (
    COLUMNS,
    default_save_path,
    run_tmprot_cli,
    score_fasta,
)


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _write_fasta(path, records):
    with open(path, "w") as f:
        for identifier, seq in records:
            f.write(f">{identifier}\n{seq}\n")


def _fake_run_writing_csv(rows, *, returncode=0, write=True):
    """subprocess.run replacement: record argv, write a synthetic tmprot CSV.

    `rows` is a list of (id, tm) written with TmProt's native header
    (Rank, ID, "Predicted Tm [°C]", Thermostable), comma-delimited.
    """
    captured = {}

    class _Proc:
        pass

    def fake_run(cmd, capture_output=True, text=True, env=None):
        captured["cmd"] = cmd
        captured["env"] = env
        in_fasta = cmd[cmd.index("-i") + 1]
        out_dir = cmd[cmd.index("-o") + 1]
        stem = os.path.splitext(os.path.basename(in_fasta))[0]
        if write and returncode == 0:
            os.makedirs(out_dir, exist_ok=True)
            with open(os.path.join(out_dir, stem + ".csv"), "w") as f:
                f.write("Rank,ID,Predicted Tm [°C],Thermostable\n")
                for rank, (identifier, tm) in enumerate(rows, 1):
                    label = "Yes" if tm > 60.0 else "No"
                    f.write(f"{rank},{identifier},{tm},{label}\n")
        p = _Proc()
        p.returncode = returncode
        p.stdout = "ok"
        p.stderr = ""
        return p

    return fake_run, captured


def test_run_tmprot_cli_command_and_output_path(monkeypatch):
    fake_run, captured = _fake_run_writing_csv([("seqA", 71.0)])
    monkeypatch.setattr(tmprot.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, [("seqA", "MKTAYIAK")])
        out_dir = os.path.join(d, "out")
        os.makedirs(out_dir)
        out_csv = run_tmprot_cli(fasta, out_dir, device="cpu")

    assert out_csv == os.path.join(out_dir, "designs.csv")
    cmd = captured["cmd"]
    assert cmd[0] == "tmprot"
    assert cmd[cmd.index("-i") + 1].endswith("designs.fasta")
    assert cmd[cmd.index("-o") + 1] == out_dir
    assert cmd[cmd.index("-d") + 1] == ","
    # device="cpu" must hide the GPU from the child process.
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == ""


def test_run_tmprot_cli_nonzero_exit_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_csv([("seqA", 60.0)], returncode=1)
    monkeypatch.setattr(tmprot.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "x.fasta")
        _write_fasta(fasta, [("seqA", "MKTAYIAK")])
        try:
            run_tmprot_cli(fasta, d)
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError on nonzero exit")


def test_run_tmprot_cli_missing_output_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_csv([("seqA", 60.0)], write=False)
    monkeypatch.setattr(tmprot.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "x.fasta")
        _write_fasta(fasta, [("seqA", "MKTAYIAK")])
        try:
            run_tmprot_cli(fasta, d)
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError when output CSV is missing")


def test_default_save_path():
    assert default_save_path("/x/y/designs.fasta") == \
        os.path.join("/x/y", "designs_tmprot.csv")
    assert default_save_path("/x/y/designs.fasta", out_suffix="tm") == \
        os.path.join("/x/y", "designs_tm.csv")


def test_score_fasta_end_to_end_with_nan_for_skipped(monkeypatch):
    """score_fasta must key by ID, keep COLUMNS order, sort by ID, name the CSV
    <stem>_tmprot.csv, drop Rank/Thermostable, and emit NaN for a FASTA ID that
    TmProt skipped (absent from its output)."""
    # TmProt returns seqA and seqC (unsorted), and SKIPS seqB (e.g. too short).
    fake_run, _ = _fake_run_writing_csv([("seqC", 48.5), ("seqA", 72.0)])
    monkeypatch.setattr(tmprot.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, [
            ("seqA", "MKTAYIAKQRQISFVK"),
            ("seqB", "MKT"),
            ("seqC", "MKTAYIAKQRQISFVK"),
        ])

        df = score_fasta(fasta)

        assert list(df.columns) == COLUMNS
        assert list(df["ID"]) == ["seqA", "seqB", "seqC"]        # sorted, full ID set
        assert os.path.isfile(os.path.join(d, "designs_tmprot.csv"))

        by_id = df.set_index("ID")
        _approx(float(by_id.loc["seqA", "tm"]), 72.0)
        _approx(float(by_id.loc["seqC", "tm"]), 48.5)
        assert np.isnan(by_id.loc["seqB", "tm"])                 # skipped -> NaN


def test_score_fasta_all_skipped_is_all_nan(monkeypatch):
    """If TmProt writes only a header (every sequence skipped), every ID is NaN."""
    fake_run, _ = _fake_run_writing_csv([])
    monkeypatch.setattr(tmprot.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, [("seqA", "MKT"), ("seqB", "MKA")])
        df = score_fasta(fasta)
    assert list(df["ID"]) == ["seqA", "seqB"]
    assert df["tm"].isna().all()


def main():
    import inspect

    class _MP:
        def __init__(self):
            self._undo = []

        def setattr(self, obj, name, value):
            self._undo.append((obj, name, getattr(obj, name)))
            setattr(obj, name, value)

        def undo(self):
            for obj, name, old in reversed(self._undo):
                setattr(obj, name, old)
            self._undo = []

    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        if "monkeypatch" in inspect.signature(t).parameters:
            mp = _MP()
            try:
                t(mp)
            finally:
                mp.undo()
        else:
            t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
