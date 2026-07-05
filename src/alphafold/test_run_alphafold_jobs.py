from __future__ import annotations

"""Self-contained tests for run_alphafold_jobs.py (AF3 job-command construction).

Run from this directory:
    cd src/alphafold && python test_run_alphafold_jobs.py
or under pytest:
    cd src/alphafold && python -m pytest test_run_alphafold_jobs.py -q

No jobs are ever submitted: we monkeypatch the module's `subprocess.run` with a stub
that RECORDS the command it was handed and returns a fake "Submitted batch job <id>"
result, then assert on the built command / --job_args string. Also tested purely:
prepare_submit_args (aurum job-name/output injection) and the duplicate-ID guard
(which raises before any submission), and skip_existing (existing .pdb -> no submit).
"""

import os
import sys
import tempfile
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pandas as pd

import alphafold.run_alphafold_jobs as raj
from alphafold.run_alphafold_jobs import prepare_submit_args, run_alphafold_jobs


class _FakeResult:
    def __init__(self, stdout):
        self.returncode = 0
        self.stdout = stdout
        self.stderr = ""


def _install_fake_subprocess(monkeypatch_calls):
    """Replace raj.subprocess.run with a recorder. Returns a restore() callable."""
    orig = raj.subprocess.run

    def fake_run(cmd, **kwargs):
        monkeypatch_calls.append(cmd)
        return _FakeResult("Submitted batch job 424242")

    raj.subprocess.run = fake_run
    return lambda: setattr(raj.subprocess, "run", orig)


def _arg_value(cmd, flag):
    """Return the token immediately following `flag` in the cmd list."""
    return cmd[cmd.index(flag) + 1]


def test_prepare_submit_args_aurum_injects_defaults():
    out = prepare_submit_args("", cluster="aurum", default_job_name="d1", working_directory="/wd")
    assert out.startswith('"') and out.endswith('"')
    inner = out[1:-1]
    assert "--job-name=AF_d1" in inner
    assert "--output=/wd/logs/%x.%j.out" in inner


def test_prepare_submit_args_respects_existing():
    out = prepare_submit_args("--job-name=MINE --output=/x.out", cluster="aurum",
                              default_job_name="d1", working_directory="/wd")
    inner = out[1:-1]
    assert "--job-name=MINE" in inner and "AF_d1" not in inner
    assert inner.count("--output=") == 1


def test_prepare_submit_args_non_aurum_passthrough():
    out = prepare_submit_args("--time=1:00:00", cluster="karolina",
                              default_job_name="d1", working_directory="/wd")
    assert out == '"--time=1:00:00"'


def test_command_and_job_args_built():
    calls = []
    restore = _install_fake_subprocess(calls)
    try:
        with tempfile.TemporaryDirectory() as wd:
            df = pd.DataFrame({"ID": ["design1"], "sequence": ["MKTAAR"]})
            run_alphafold_jobs(
                df=df, working_directory=wd, cluster="aurum",
                protein_column_names=[("ID", "sequence")],
                skip_existing=True,
            )
        assert len(calls) == 1, calls
        cmd = calls[0]
        assert cmd[0] == "bash"
        assert _arg_value(cmd, "--cluster") == "aurum"
        assert _arg_value(cmd, "--job_name") == "alphafold"
        job_args = _arg_value(cmd, "--job_args")
        assert f"--working_directory {wd}" in job_args
        assert "--sequence_id design1" in job_args
        assert "--proteins design1 MKTAAR" in job_args
        assert "--model_seeds 42" in job_args
        # apo (no ligands/ions) -> those flags absent
        assert "--ligands" not in job_args and "--ions" not in job_args
    finally:
        restore()


def test_ligands_and_ions_in_job_args():
    calls = []
    restore = _install_fake_subprocess(calls)
    try:
        with tempfile.TemporaryDirectory() as wd:
            df = pd.DataFrame({
                "ID": ["d1"], "sequence": ["MK"],
                "lig_id": ["LIG"], "lig_smi": ["CCO"],
                "ion_id": ["MG1"], "ion_ccd": ["MG"],
            })
            run_alphafold_jobs(
                df=df, working_directory=wd, cluster="aurum",
                protein_column_names=[("ID", "sequence")],
                ligand_column_names=[("lig_id", "lig_smi")],
                ion_column_names=[("ion_id", "ion_ccd")],
            )
        job_args = _arg_value(calls[0], "--job_args")
        assert "--ligands LIG CCO" in job_args
        assert "--ions MG1 MG" in job_args
    finally:
        restore()


def test_duplicate_ids_raise():
    calls = []
    restore = _install_fake_subprocess(calls)
    try:
        with tempfile.TemporaryDirectory() as wd:
            df = pd.DataFrame({"ID": ["dup", "dup"], "sequence": ["MK", "GG"]})
            try:
                run_alphafold_jobs(
                    df=df, working_directory=wd, cluster="aurum",
                    protein_column_names=[("ID", "sequence")],
                )
            except ValueError:
                assert calls == [], "must raise before submitting anything"
                return
            raise AssertionError("expected ValueError on duplicate folding IDs")
    finally:
        restore()


def test_skip_existing_pdb_not_submitted():
    calls = []
    restore = _install_fake_subprocess(calls)
    try:
        with tempfile.TemporaryDirectory() as wd:
            os.makedirs(os.path.join(wd, "structs"))
            # design "have" already has a pdb; "need" does not.
            Path(os.path.join(wd, "structs", "have.pdb")).write_text("X")
            df = pd.DataFrame({"ID": ["have", "need"], "sequence": ["MK", "GG"]})
            run_alphafold_jobs(
                df=df, working_directory=wd, cluster="aurum",
                protein_column_names=[("ID", "sequence")],
                skip_existing=True,
            )
        assert len(calls) == 1, calls
        assert "--sequence_id need" in _arg_value(calls[0], "--job_args")
    finally:
        restore()


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
