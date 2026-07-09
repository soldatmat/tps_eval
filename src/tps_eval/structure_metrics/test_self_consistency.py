from __future__ import annotations

"""Self-contained tests for self_consistency.py (scRMSD — PURE-PYTHON parts).

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_self_consistency.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_self_consistency.py -q

Neither ProteinMPNN (torch) nor ESMFold is invoked. ESMFold is dependency-injected as
`fold_fn`, so self_consistency_for_structure is driven with a FAKE folder; ProteinMPNN's
sampling is exercised by monkeypatching subprocess.run to drop a synthetic ProteinMPNN
`.fa` (so we test the argv AND the native-drop / '/'-join parse). The Cα-RMSD math is
checked against closed-form synthetic PDBs:
  * identical / translated / rotated backbones -> RMSD 0 (superposition removes rigid moves),
  * a symmetric 4-atom set scaled x2 -> optimal rotation is identity, RMSD == 1.0 exactly,
  * residue-count mismatch -> aligns the leading min(len) atoms,
and the orchestration (n_samples counting, min/mean, NaN-on-all-failure), collection,
and default CSV naming.

self_consistency_dir itself is NOT run end-to-end (it imports torch-backed esmfold at
call time). That single path is NEEDS-AURUM; everything it composes is covered here.
"""

import os
import tempfile

import numpy as np
import pandas as pd

import tps_eval.structure_metrics.self_consistency as self_consistency
from tps_eval.structure_metrics.self_consistency import (
    COLUMNS,
    _ca_atoms,
    _ca_rmsd,
    _chain_ids,
    _collect_structures,
    _default_save_path,
    _sample_sequences,
    self_consistency_for_structure,
)


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _write_ca_pdb(path, coords, chain="A", resname="ALA"):
    """Minimal PDB, one CA per residue at the given coords."""
    with open(path, "w") as fh:
        for i, (x, y, z) in enumerate(coords, start=1):
            fh.write(
                f"ATOM  {i:>5d}  CA  {resname} {chain}{i:>4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        fh.write("END\n")


def test_ca_atoms_order_and_count():
    coords = [(0, 0, 0), (1, 0, 0), (2, 0, 0)]
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "s.pdb")
        _write_ca_pdb(p, coords)
        cas = _ca_atoms(p)
        assert len(cas) == 3
        got = np.array([a.get_coord() for a in cas])
        np.testing.assert_allclose(got, np.array(coords, float), atol=1e-3)


def test_chain_ids_only_polymer_chains():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "s.pdb")
        _write_ca_pdb(p, [(0, 0, 0), (3.8, 0, 0)], chain="A")
        assert _chain_ids(p) == ["A"]


def test_ca_rmsd_identical_zero():
    coords = [(0, 0, 0), (3.8, 0, 0), (7.6, 0, 0), (11.4, 0, 0)]
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.pdb")
        b = os.path.join(d, "b.pdb")
        _write_ca_pdb(a, coords)
        _write_ca_pdb(b, coords)
        _approx(_ca_rmsd(a, b), 0.0)


def test_ca_rmsd_rigid_move_zero():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(12, 3)) * 5.0
    # A rigid rotation about z by 37 deg + a translation must superpose to ~0 RMSD.
    th = np.deg2rad(37.0)
    R = np.array([[np.cos(th), -np.sin(th), 0],
                  [np.sin(th), np.cos(th), 0],
                  [0, 0, 1.0]])
    moved = coords @ R.T + np.array([100.0, -20.0, 5.0])
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.pdb")
        b = os.path.join(d, "b.pdb")
        _write_ca_pdb(a, coords)
        _write_ca_pdb(b, moved)
        _approx(_ca_rmsd(a, b), 0.0, tol=1e-3)


def test_ca_rmsd_scaled_symmetric_set_is_one():
    """Closed form: for the symmetric unit set {(±1,0,0),(0,±1,0)} the optimal
    rotation onto 2× itself is the identity, so residual == the original vectors
    and RMSD == sqrt(mean|v|^2) == 1.0 exactly."""
    base = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], float)
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.pdb")
        b = os.path.join(d, "b.pdb")
        _write_ca_pdb(a, base)
        _write_ca_pdb(b, base * 2.0)
        _approx(_ca_rmsd(a, b), 1.0, tol=1e-4)


def test_ca_rmsd_length_mismatch_aligns_leading():
    coords = [(0, 0, 0), (3.8, 0, 0), (7.6, 0, 0)]
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.pdb")
        b = os.path.join(d, "b.pdb")
        _write_ca_pdb(a, coords)
        _write_ca_pdb(b, coords + [(11.4, 0, 0), (15.2, 0, 0)])  # 2 extra residues
        _approx(_ca_rmsd(a, b), 0.0, tol=1e-3)     # leading 3 align exactly


def test_ca_rmsd_no_polymer_ca_is_nan():
    """A structure that parses but has no polymer Cα (e.g. HETATM-only) hits the
    `if not ref or not mob` guard and yields NaN rather than a bogus RMSD."""
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.pdb")
        b = os.path.join(d, "b.pdb")
        _write_ca_pdb(a, [(0, 0, 0)])
        with open(b, "w") as fh:  # a lone HETATM: model parses, but no polymer CA
            fh.write(
                "HETATM    1 MG    MG A 900       0.000   0.000   0.000  1.00  0.00          MG\n"
            )
            fh.write("END\n")
        assert np.isnan(_ca_rmsd(a, b))


def test_collect_structures_flat_and_af3():
    with tempfile.TemporaryDirectory() as d:
        for name in ("a.pdb", "a.cif", "b.cif"):
            open(os.path.join(d, name), "w").close()
        structures, mode = _collect_structures(d)
        assert mode == "flat"
        assert list(structures.keys()) == ["a", "b"]
        assert structures["a"].endswith("a.pdb")
    with tempfile.TemporaryDirectory() as d:
        job = os.path.join(d, "seq1")
        os.makedirs(job)
        open(os.path.join(job, "seq1_model.cif"), "w").close()
        structures, mode = _collect_structures(d)
        assert mode == "af3"


def test_default_save_path():
    assert _default_save_path("/x/y/structs") == \
        os.path.join("/x/y", "structs_self_consistency.csv")


def test_sample_sequences_command_and_parse(monkeypatch):
    """ProteinMPNN sampling: assert the argv, and that the native (first) record is
    dropped, the N sampled records kept, and multi-chain '/' separators removed."""
    captured = {}

    class _Proc:
        returncode = 0
        stdout = "ok"
        stderr = ""

    def fake_run(cmd, capture_output=True, text=True):
        captured["cmd"] = cmd
        out_folder = cmd[cmd.index("--out_folder") + 1]
        pdb_path = cmd[cmd.index("--pdb_path") + 1]
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        seqs_dir = os.path.join(out_folder, "seqs")
        os.makedirs(seqs_dir, exist_ok=True)
        with open(os.path.join(seqs_dir, stem + ".fa"), "w") as fh:
            fh.write(">native, T=0.1, score=1.0\nMKLNATIVE\n")     # native -> dropped
            fh.write(">sample_1, T=0.1\nAAAA\n")
            fh.write(">sample_2, T=0.1\nCCC/DDD\n")                  # multichain
        return _Proc()

    monkeypatch.setattr(self_consistency.subprocess, "run", fake_run)
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "mpnnout")
        os.makedirs(out)
        got = _sample_sequences(
            os.path.join(d, "design3.pdb"), out,
            num_seqs=2, sampling_temp=0.1, model_name="v_48_020", seed=0,
            pdb_path_chains="A",
        )
    assert got == ["AAAA", "CCCDDD"]     # native dropped; '/' removed

    cmd = captured["cmd"]
    assert str(self_consistency.PROTEINMPNN_DIR / "protein_mpnn_run.py") in cmd
    assert cmd[cmd.index("--num_seq_per_target") + 1] == "2"
    assert cmd[cmd.index("--sampling_temp") + 1] == "0.1"
    assert cmd[cmd.index("--model_name") + 1] == "v_48_020"
    assert cmd[cmd.index("--pdb_path_chains") + 1] == "A"
    assert cmd[cmd.index("--batch_size") + 1] == "1"


def test_self_consistency_for_structure_fake_fold(monkeypatch):
    """Drive the sample->refold->RMSD loop with a fake folder that returns the ref
    backbone verbatim (=> RMSD 0) for 2 of 3 samples and raises for the third."""
    coords = [(0, 0, 0), (3.8, 0, 0), (7.6, 0, 0), (11.4, 0, 0)]
    with tempfile.TemporaryDirectory() as d:
        ref = os.path.join(d, "designA.pdb")
        _write_ca_pdb(ref, coords)
        ref_text = open(ref).read()

        monkeypatch.setattr(self_consistency, "_sample_sequences",
                            lambda *a, **k: ["S1", "S2", "S3"])

        calls = {"n": 0}

        def fake_fold(seq):
            calls["n"] += 1
            if seq == "S3":
                raise RuntimeError("fold failed")
            return ref_text          # identical backbone -> RMSD 0

        workdir = os.path.join(d, "work")
        os.makedirs(workdir)
        res = self_consistency_for_structure(
            ref, fold_fn=fake_fold, num_seqs=3, sampling_temp=0.1,
            model_name="v_48_020", seed=0, workdir=workdir,
        )
    _approx(res["sc_rmsd_min"], 0.0, tol=1e-3)
    _approx(res["sc_rmsd_mean"], 0.0, tol=1e-3)
    assert res["n_samples"] == 2      # S3 raised -> not counted


def test_self_consistency_for_structure_no_chains_nan():
    """A structure with no polymer CA atoms yields a graceful NaN row."""
    with tempfile.TemporaryDirectory() as d:
        empty = os.path.join(d, "empty.pdb")
        with open(empty, "w") as fh:  # HETATM only: parses to a model with no polymer chain
            fh.write(
                "HETATM    1 MG    MG A 900       0.000   0.000   0.000  1.00  0.00          MG\n"
            )
            fh.write("END\n")
        workdir = os.path.join(d, "work")
        os.makedirs(workdir)
        res = self_consistency_for_structure(
            empty, fold_fn=lambda s: "", num_seqs=2, sampling_temp=0.1,
            model_name="v_48_020", seed=0, workdir=workdir,
        )
    assert np.isnan(res["sc_rmsd_min"])
    assert res["n_samples"] == 0


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
