from __future__ import annotations

"""Self-contained tests for proteinmpnn_score.py (ProteinMPNN NLL — PURE-PYTHON parts).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_proteinmpnn_score.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_proteinmpnn_score.py -q

ProteinMPNN (torch) is NEVER invoked. We monkeypatch subprocess.run to (a) capture the
command that WOULD be launched and (b) drop a synthetic score_only/<stem>_pdb.npz, then
assert:
  * score_pdb builds the right argv (--score_only 1, --model_name, --pdb_path/--out_folder,
    --batch_size 1, --num_seq_per_target, --seed, --backbone_noise) and averages the npz
    `global_score` / `score` arrays into proteinmpnn_nll / proteinmpnn_score_designed,
  * a nonzero exit and a missing npz both raise RuntimeError,
  * flat/af3 structure collection + ID keying + default CSV naming,
  * score_dir DataFrame assembly, column order, sort, NaN-on-failure, CSV naming.
"""

import os
import tempfile

import numpy as np
import pandas as pd

import proteinmpnn_score
from proteinmpnn_score import (
    COLUMNS,
    _collect_structures,
    _default_save_path,
    score_dir,
    score_pdb,
)


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _fake_run_writing_npz(*, global_score, score, returncode=0, write=True):
    """Build a subprocess.run replacement that records the argv and writes an npz."""
    captured = {}

    class _Proc:
        pass

    def fake_run(cmd, capture_output=True, text=True):
        captured["cmd"] = cmd
        out_folder = cmd[cmd.index("--out_folder") + 1]
        pdb_path = cmd[cmd.index("--pdb_path") + 1]
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        if write and returncode == 0:
            score_dir_ = os.path.join(out_folder, "score_only")
            os.makedirs(score_dir_, exist_ok=True)
            np.savez(os.path.join(score_dir_, stem + "_pdb.npz"),
                     global_score=np.asarray(global_score, float),
                     score=np.asarray(score, float))
        p = _Proc()
        p.returncode = returncode
        p.stdout = "ok"
        p.stderr = ""
        return p

    return fake_run, captured


def test_score_pdb_command_and_npz_parsing(monkeypatch):
    fake_run, captured = _fake_run_writing_npz(
        global_score=[1.0, 2.0, 3.0], score=[4.0, 6.0])
    monkeypatch.setattr(proteinmpnn_score.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "out")
        os.makedirs(out)
        res = score_pdb(os.path.join(d, "design7.pdb"), out,
                        model_name="v_48_030", seed=3, num_passes=1, backbone_noise=0.0)
    _approx(res["proteinmpnn_nll"], 2.0)                # mean(global_score)
    _approx(res["proteinmpnn_score_designed"], 5.0)     # mean(score)

    cmd = captured["cmd"]
    assert cmd[0].endswith("python") or os.path.basename(cmd[0]).startswith("python") \
        or cmd[0] == proteinmpnn_score.sys.executable
    assert str(proteinmpnn_score.PROTEINMPNN_DIR / "protein_mpnn_run.py") in cmd
    assert cmd[cmd.index("--score_only") + 1] == "1"
    assert cmd[cmd.index("--pdb_path") + 1].endswith("design7.pdb")
    assert cmd[cmd.index("--batch_size") + 1] == "1"
    assert cmd[cmd.index("--seed") + 1] == "3"
    assert cmd[cmd.index("--num_seq_per_target") + 1] == "1"
    # Regression: --model_name must be forwarded (was silently dropped before).
    assert cmd[cmd.index("--model_name") + 1] == "v_48_030"


def test_score_pdb_nonzero_exit_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_npz(global_score=[1.0], score=[1.0], returncode=1)
    monkeypatch.setattr(proteinmpnn_score.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        try:
            score_pdb(os.path.join(d, "x.pdb"), d)
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError on nonzero exit")


def test_score_pdb_missing_npz_raises(monkeypatch):
    fake_run, _ = _fake_run_writing_npz(global_score=[1.0], score=[1.0], write=False)
    monkeypatch.setattr(proteinmpnn_score.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        try:
            score_pdb(os.path.join(d, "x.pdb"), d)
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError when npz is missing")


def test_collect_structures_flat_pdb_wins():
    with tempfile.TemporaryDirectory() as d:
        for name in ("a.pdb", "a.cif", "b.cif"):
            open(os.path.join(d, name), "w").close()
        structures, mode = _collect_structures(d)
        assert mode == "flat"
        assert list(structures.keys()) == ["a", "b"]     # OrderedDict, sorted
        assert structures["a"].endswith("a.pdb")
        assert structures["b"].endswith("b.cif")


def test_collect_structures_af3_layout():
    with tempfile.TemporaryDirectory() as d:
        job = os.path.join(d, "seqX")
        os.makedirs(job)
        open(os.path.join(job, "seqX_model.cif"), "w").close()
        structures, mode = _collect_structures(d)
        assert mode == "af3"
        assert structures["seqX"].endswith(os.path.join("seqX", "seqX_model.cif"))


def test_collect_structures_missing_dir_is_empty():
    structures, mode = _collect_structures("/no/such/dir/anywhere")
    assert len(structures) == 0
    assert mode == "flat"


def test_default_save_path():
    assert _default_save_path("/x/y/structs") == \
        os.path.join("/x/y", "structs_proteinmpnn_score.csv")
    assert _default_save_path("/x/y/structs/") == \
        os.path.join("/x/y", "structs_proteinmpnn_score.csv")


def test_score_dir_end_to_end_with_nan_on_failure(monkeypatch):
    """score_dir must key by ID, keep COLUMNS order, sort by ID, name the CSV
    <structs_dir>_proteinmpnn_score.csv, and emit a NaN row when one structure fails."""
    def fake_score_pdb(path, out_folder, *, model_name, seed, backbone_noise):
        stem = os.path.splitext(os.path.basename(path))[0]
        if stem == "bad":
            raise RuntimeError("mpnn failed")
        return {"proteinmpnn_nll": 1.25 if stem == "good" else 2.5,
                "proteinmpnn_score_designed": 0.5}

    monkeypatch.setattr(proteinmpnn_score, "score_pdb", fake_score_pdb)
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        for name in ("good.pdb", "bad.pdb", "aaa.pdb"):
            open(os.path.join(structs, name), "w").close()

        df = score_dir(structs)

        assert list(df.columns) == COLUMNS
        assert list(df["ID"]) == ["aaa", "bad", "good"]        # sorted
        assert os.path.isfile(structs + "_proteinmpnn_score.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["proteinmpnn_nll"]), 1.25)

        bad = df.set_index("ID").loc["bad"]
        assert np.isnan(bad["proteinmpnn_nll"])
        assert np.isnan(bad["proteinmpnn_score_designed"])


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
