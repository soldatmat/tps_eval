# -*- coding: utf-8 -*-
"""Self-contained tests for aggregation.py (Aggrescan3D wrapper — PURE-PYTHON parts).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_aggregation.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_aggregation.py -q

The real Aggrescan3D binary (Python-2.7 vendored tool) is NEVER invoked. We test:
  * A3D.csv parsing (_parse_a3d_csv) on synthetic fixtures (header, junk rows, empty),
  * the per-ID scalar reduction (_summarize) with closed-form arrays incl. the
    positive-only sum and the empty -> NaN contract,
  * the command that _run_a3d_static WOULD run (subprocess.Popen monkeypatched: assert
    static mode, i.e. NO --dynamic flag, plus the -i/-w args) and its failure handling,
  * the flat/af3 structure collection + ID keying + default CSV naming,
  * end-to-end extract_aggregation_dir DataFrame assembly, ID keying, CSV naming and the
    NaN-on-failure contract (A3D itself stubbed to write a synthetic A3D.csv / to raise).

NOTE: aggregation.py must stay Python-2.7-compatible (it runs in the aggrescan3d env),
so this test avoids py3-only constructs in anything it feeds back into that module.
"""

import os
import tempfile

import numpy as np
import pandas as pd

import tps_eval.structure_metrics.aggregation as aggregation
from tps_eval.structure_metrics.aggregation import (
    COLUMNS,
    _collect_structures,
    _default_save_path,
    _parse_a3d_csv,
    _run_a3d_static,
    _summarize,
    extract_aggregation_dir,
)


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, "%s != %s (tol %s)" % (a, b, tol)


def _write_a3d_csv(path, scores):
    """Write a synthetic A3D.csv (protein,chain,residue,residue_name,score)."""
    with open(path, "w") as fh:
        fh.write("protein,chain,residue,residue_name,score\n")
        for i, s in enumerate(scores, start=1):
            fh.write("prot,A,%d,ALA,%s\n" % (i, s))


def test_parse_a3d_csv_basic():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "A3D.csv")
        _write_a3d_csv(p, [1.5, -0.5, 0.0, 2.25])
        got = _parse_a3d_csv(p)
        assert got.dtype == float
        np.testing.assert_allclose(got, [1.5, -0.5, 0.0, 2.25])


def test_parse_a3d_csv_skips_junk_and_short_rows():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "A3D.csv")
        with open(p, "w") as fh:
            fh.write("protein,chain,residue,residue_name,score\n")
            fh.write("prot,A,1,ALA,1.0\n")
            fh.write("too,short,row\n")             # <5 cols -> skipped
            fh.write("prot,A,2,ALA,not_a_float\n")   # unparsable score -> skipped
            fh.write("prot,A,3,ALA,-2.0\n")
        got = _parse_a3d_csv(p)
        np.testing.assert_allclose(got, [1.0, -2.0])


def test_parse_a3d_csv_header_only_is_empty():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "A3D.csv")
        _write_a3d_csv(p, [])
        got = _parse_a3d_csv(p)
        assert got.shape == (0,)


def test_summarize_closed_form():
    m = _summarize(np.array([2.0, -1.0, 3.0, -4.0]))
    _approx(m["a3d_total_score"], 0.0)
    _approx(m["a3d_avg_score"], 0.0)
    _approx(m["a3d_max_score"], 3.0)
    _approx(m["a3d_min_score"], -4.0)
    _approx(m["a3d_total_pos_score"], 5.0)   # 2 + 3, negatives excluded
    assert m["n_residues"] == 4


def test_summarize_all_negative_pos_score_zero():
    m = _summarize(np.array([-1.0, -2.0]))
    _approx(m["a3d_total_pos_score"], 0.0)   # empty positive set sums to 0.0
    _approx(m["a3d_max_score"], -1.0)


def test_summarize_empty_is_nan():
    m = _summarize(np.empty((0,)))
    assert m["n_residues"] == 0
    for k in ("a3d_avg_score", "a3d_total_score", "a3d_max_score",
              "a3d_min_score", "a3d_total_pos_score"):
        assert np.isnan(m[k]), k


def test_run_a3d_static_builds_static_command(monkeypatch):
    """The built command must be STATIC (no --dynamic) and reference the given
    pdb/work paths; on success it returns the produced A3D.csv path."""
    captured = {}

    class _FakeProc(object):
        returncode = 0

        def communicate(self):
            return (b"aggrescan done\n", None)

    def fake_popen(cmd, stdout=None, stderr=None):
        captured["cmd"] = cmd
        # A3D writes A3D.csv into -w work_dir; emulate that.
        w = cmd[cmd.index("-w") + 1]
        if not os.path.isdir(w):
            os.makedirs(w)
        _write_a3d_csv(os.path.join(w, "A3D.csv"), [1.0, 2.0])
        return _FakeProc()

    monkeypatch.setattr(aggregation.subprocess, "Popen", fake_popen)
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, "run")
        csv_path = _run_a3d_static(os.path.join(d, "x.pdb"), work)
        assert csv_path == os.path.join(work, "A3D.csv")
        assert os.path.isfile(csv_path)

    cmd = captured["cmd"]
    assert cmd[0] == "aggrescan"
    assert "--dynamic" not in cmd and "-d" not in cmd     # STATIC mode only
    assert cmd[cmd.index("-i") + 1].endswith("x.pdb")
    assert cmd[cmd.index("-w") + 1] == work
    assert "--overwrite" in cmd


def test_run_a3d_static_failure_raises(monkeypatch):
    class _FailProc(object):
        returncode = 1

        def communicate(self):
            return (b"boom\n", None)

    monkeypatch.setattr(aggregation.subprocess, "Popen",
                        lambda *a, **k: _FailProc())
    with tempfile.TemporaryDirectory() as d:
        try:
            _run_a3d_static(os.path.join(d, "x.pdb"), os.path.join(d, "run"))
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected RuntimeError on nonzero exit / missing A3D.csv")


def test_collect_structures_flat_pdb_wins(monkeypatch):
    with tempfile.TemporaryDirectory() as d:
        for name in ("a.pdb", "a.cif", "b.cif"):
            open(os.path.join(d, name), "w").close()
        pairs, mode = _collect_structures(d)
        assert mode == "flat"
        m = dict(pairs)
        assert set(m) == {"a", "b"}
        assert m["a"].endswith("a.pdb")   # .pdb wins over .cif for the same stem
        assert m["b"].endswith("b.cif")
        # sorted by ID
        assert [k for k, _ in pairs] == ["a", "b"]


def test_collect_structures_af3_layout():
    with tempfile.TemporaryDirectory() as d:
        job = os.path.join(d, "seq1")
        os.makedirs(job)
        open(os.path.join(job, "seq1_model.cif"), "w").close()
        pairs, mode = _collect_structures(d)
        assert mode == "af3"
        assert dict(pairs)["seq1"].endswith(os.path.join("seq1", "seq1_model.cif"))


def test_default_save_path():
    p = _default_save_path("/tmp/foo/my_structs")
    assert p == os.path.join("/tmp/foo", "my_structs_aggregation.csv")
    # trailing sep is stripped
    p2 = _default_save_path("/tmp/foo/my_structs/")
    assert p2 == os.path.join("/tmp/foo", "my_structs_aggregation.csv")


def test_extract_dir_end_to_end_with_nan_on_failure(monkeypatch):
    """End-to-end orchestration with A3D stubbed: good structures score, one
    structure fails -> NaN row (processing continues); ID keying, column order,
    sort and default CSV naming are checked."""
    # Stub the two A3D-touching steps so no external binary / Biopython PDB is needed.
    monkeypatch.setattr(aggregation, "_prepare_pdb",
                        lambda src, out: (open(out, "w").close() or out))

    def fake_run(pdb_path, work_dir):
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        if stem == "bad":
            raise RuntimeError("A3D blew up")
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)
        csv_path = os.path.join(work_dir, "A3D.csv")
        _write_a3d_csv(csv_path, [1.0, -1.0, 3.0] if stem == "good" else [0.5])
        return csv_path

    monkeypatch.setattr(aggregation, "_run_a3d_static", fake_run)

    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        for name in ("good.pdb", "bad.pdb", "z_other.pdb"):
            open(os.path.join(structs, name), "w").close()

        df = extract_aggregation_dir(structs)

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "bad", "z_other"}
        # sorted by ID
        assert list(df["ID"]) == ["bad", "good", "z_other"]
        assert os.path.isfile(structs + "_aggregation.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["a3d_total_score"]), 3.0)
        _approx(float(good["a3d_total_pos_score"]), 4.0)   # 1 + 3
        assert int(good["n_residues"]) == 3

        bad = df.set_index("ID").loc["bad"]
        assert np.isnan(bad["a3d_avg_score"])
        assert int(bad["n_residues"]) == 0

        other = df.set_index("ID").loc["z_other"]
        _approx(float(other["a3d_total_score"]), 0.5)
        assert int(other["n_residues"]) == 1


def main():
    import inspect
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    # Minimal monkeypatch shim so the file runs standalone (no pytest needed).
    class _MP(object):
        def __init__(self):
            self._undo = []

        def setattr(self, obj, name, value):
            old = getattr(obj, name)
            self._undo.append((obj, name, old))
            setattr(obj, name, value)

        def undo(self):
            for obj, name, old in reversed(self._undo):
                setattr(obj, name, old)
            self._undo = []

    for t in tests:
        needs_mp = "monkeypatch" in inspect.signature(t).parameters
        if needs_mp:
            mp = _MP()
            try:
                t(mp)
            finally:
                mp.undo()
        else:
            t()
        print("  ok  %s" % t.__name__)
    print("All %d tests passed." % len(tests))


if __name__ == "__main__":
    main()
