from __future__ import annotations

"""Self-contained tests for catapro.py (CataPro kinetics — PURE-PYTHON parts).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/sequence_metrics && python test_catapro.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_catapro.py -q

CataPro (torch/ProtT5/MolT5) is NEVER invoked. We monkeypatch subprocess.run to
(a) capture the argv + cwd that WOULD be launched and (b) drop a synthetic native
CataPro output CSV, then assert:
  * score_fasta builds the predict.py argv (-inp_fpath/-out_fpath/-model_dpath/
    -device) and runs it with cwd = the CataPro inference dir,
  * the native log10 columns are exponentiated to absolute catapro_kcat/km/kcat_km,
    IDs are recovered by stripping the "_wild" tag, keyed by ID, COLUMNS order, sorted,
    and written to <fasta_stem>_<out_suffix>.csv,
  * a sequence CataPro drops is reindexed back in as a NaN row,
  * an unknown substrate (no SMILES) yields all-NaN rows WITHOUT invoking CataPro,
  * a nonzero exit and a missing output file both raise RuntimeError,
  * resolve_smiles + log10_to_absolute helpers.
"""

import os
import tempfile

import numpy as np
import pandas as pd

import catapro
from catapro import (
    COLUMNS,
    NATIVE_ID,
    NATIVE_KCAT,
    NATIVE_KM,
    NATIVE_KCAT_KM,
    log10_to_absolute,
    resolve_smiles,
    score_fasta,
)

# Fixed log10 predictions -> absolute: 10**2=100, 10**0=1, 10**-1=0.1
_KCAT_LOG10 = 2.0
_KM_LOG10 = 0.0
_KCAT_KM_LOG10 = -1.0


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _write_fasta(path, ids):
    with open(path, "w") as f:
        for i, seq_id in enumerate(ids):
            f.write(f">{seq_id} description_{i}\nMKAILVTDPRSTQWACDEFGHIKLMNPQ\n")


def _fake_run_writing_native(*, returncode=0, write=True, skip_ids=()):
    """subprocess.run replacement: record (cmd, cwd), write a synthetic native CSV."""
    captured = {}

    class _Proc:
        pass

    def fake_run(cmd, cwd=None, capture_output=True, text=True):
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        inp = cmd[cmd.index("-inp_fpath") + 1]
        out = cmd[cmd.index("-out_fpath") + 1]
        if write and returncode == 0:
            inp_df = pd.read_csv(inp, index_col=0)
            rows = []
            for enzyme_id, smiles in zip(inp_df["Enzyme_id"], inp_df["smiles"]):
                if enzyme_id in skip_ids:
                    continue
                rows.append({
                    NATIVE_ID: f"{enzyme_id}_wild",
                    "smiles": smiles,
                    NATIVE_KCAT: _KCAT_LOG10,
                    NATIVE_KM: _KM_LOG10,
                    NATIVE_KCAT_KM: _KCAT_KM_LOG10,
                })
            pd.DataFrame(rows).to_csv(out)  # native predict.py writes an index col
        p = _Proc()
        p.returncode = returncode
        p.stdout = "ok"
        p.stderr = ""
        return p

    return fake_run, captured


def test_resolve_smiles_and_conversion():
    assert resolve_smiles("FPP") == catapro.SUBSTRATE_SMILES["FPP"]
    assert resolve_smiles("fpp") == catapro.SUBSTRATE_SMILES["FPP"]      # case-insensitive
    assert resolve_smiles("NOPE") is None                               # unknown -> None
    assert resolve_smiles("NOPE", smiles="CCO") == "CCO"                # explicit wins
    assert resolve_smiles(None) is None
    _approx(log10_to_absolute(2.0), 100.0)
    _approx(log10_to_absolute(0.0), 1.0)
    _approx(log10_to_absolute(-1.0), 0.1)


def test_score_fasta_reshape_and_command(monkeypatch):
    fake_run, captured = _fake_run_writing_native()
    monkeypatch.setattr(catapro.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, ["d2", "d1", "d3"])
        df = score_fasta(fasta, "FPP", device="cuda")
        assert os.path.isfile(os.path.join(d, "designs_catapro.csv"))

    assert list(df.columns) == COLUMNS
    assert list(df["ID"]) == ["d1", "d2", "d3"]            # sorted
    row = df.set_index("ID").loc["d1"]
    _approx(float(row["catapro_kcat"]), 100.0)
    _approx(float(row["catapro_km"]), 1.0)
    _approx(float(row["catapro_kcat_km"]), 0.1)
    assert (df["catapro_substrate"] == "FPP").all()

    cmd = captured["cmd"]
    assert str(catapro.CATAPRO_PREDICT) in cmd
    assert cmd[cmd.index("-device") + 1] == "cuda"
    assert cmd[cmd.index("-inp_fpath") + 1].endswith(".csv")
    assert cmd[cmd.index("-out_fpath") + 1].endswith(".csv")
    # predict.py uses flat imports -> must run with cwd = the inference dir.
    assert os.path.abspath(captured["cwd"]) == os.path.abspath(str(catapro.CATAPRO_INFERENCE_DIR))


def test_out_suffix_names_the_file(monkeypatch):
    fake_run, _ = _fake_run_writing_native()
    monkeypatch.setattr(catapro.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "marts.fasta")
        _write_fasta(fasta, ["m1"])
        score_fasta(fasta, "GGPP", out_suffix="catapro_GGPP")
        assert os.path.isfile(os.path.join(d, "marts_catapro_GGPP.csv"))


def test_dropped_sequence_becomes_nan_row(monkeypatch):
    fake_run, _ = _fake_run_writing_native(skip_ids=("d2",))
    monkeypatch.setattr(catapro.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, ["d1", "d2", "d3"])
        df = score_fasta(fasta, "FPP")
    assert list(df["ID"]) == ["d1", "d2", "d3"]
    dropped = df.set_index("ID").loc["d2"]
    assert np.isnan(dropped["catapro_kcat"])
    assert np.isnan(dropped["catapro_km"])
    assert np.isnan(dropped["catapro_kcat_km"])
    # substrate label is still recorded on the reindexed NaN row
    assert dropped["catapro_substrate"] == "FPP"


def test_unknown_substrate_all_nan_without_invoking(monkeypatch):
    def explode(*a, **k):
        raise AssertionError("CataPro must NOT be invoked when no SMILES is known")
    monkeypatch.setattr(catapro.subprocess, "run", explode)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, ["d1", "d2"])
        df = score_fasta(fasta, "2xGGPP")   # no SMILES in the map
        assert os.path.isfile(os.path.join(d, "designs_catapro.csv"))
    assert list(df.columns) == COLUMNS
    assert list(df["ID"]) == ["d1", "d2"]
    assert df["catapro_kcat"].isna().all()
    assert df["catapro_km"].isna().all()
    assert df["catapro_kcat_km"].isna().all()
    assert (df["catapro_substrate"] == "2XGGPP").all()


def test_nonzero_exit_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_native(returncode=1)
    monkeypatch.setattr(catapro.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, ["d1"])
        try:
            score_fasta(fasta, "FPP")
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError on nonzero exit")


def test_missing_output_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_native(write=False)
    monkeypatch.setattr(catapro.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "designs.fasta")
        _write_fasta(fasta, ["d1"])
        try:
            score_fasta(fasta, "FPP")
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError when output file is missing")


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
