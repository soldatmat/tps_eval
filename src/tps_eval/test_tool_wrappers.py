"""Repo-integrity checks on the shell wrappers under ``scripts/``.

Run: python -m pytest src/tps_eval/test_tool_wrappers.py

Every eval tool activates its OWN conda env (``paths.sh``: aggrescan3d, pocket,
tmprot, catapro, enzyme_explorer, esmfold, ...) and then runs
``python -m tps_eval.<subdir>.run_<tool>``. Only the main ``tps_eval`` env ever
gets ``pip install -e .``, so the wrappers must put the in-repo package on
PYTHONPATH themselves — otherwise a repo pull silently leaves the satellite envs
importing a stale (or absent) copy and five tools die with
``ModuleNotFoundError: tps_eval`` on the next submission round.

These are static checks: they guard the invariant for wrappers added later.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WRAPPER_DIR = REPO_ROOT / "scripts" / "tool_wrappers"

# The exact line the wrappers use to expose the in-repo package.
PYTHONPATH_LINE = 'export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"'

# Top-level scripts that also activate an env and run the package.
TOP_LEVEL_SCRIPTS = [
    "run_selection.sh",
    "run_eval_pipeline_continuation.sh",
    "run_prepare_order.sh",
    "run_alphafold_fanout.sh",
]

_CD_TO_ROOT = re.compile(r'^cd "\$SCRIPT_DIR/(\.\./)?\.\."( \|\| exit 1)?$')


def _scripts_running_the_package():
    paths = sorted(WRAPPER_DIR.glob("run_*.sh"))
    paths += [REPO_ROOT / "scripts" / name for name in TOP_LEVEL_SCRIPTS]
    return [p for p in paths if "python -m tps_eval" in p.read_text()]


def test_wrapper_dir_is_where_we_think_it_is():
    assert WRAPPER_DIR.is_dir(), WRAPPER_DIR
    assert len(_scripts_running_the_package()) > 20


def test_every_wrapper_exposes_the_in_repo_package():
    missing = [p.name for p in _scripts_running_the_package() if PYTHONPATH_LINE not in p.read_text()]
    assert not missing, (
        "these scripts run `python -m tps_eval` without putting the repo's src/ on "
        "PYTHONPATH, so they break in any env lacking an up-to-date editable install: "
        + ", ".join(missing)
    )


def test_pythonpath_export_comes_after_cd_to_repo_root_and_before_python():
    """`$(pwd)` must resolve to the repo root, and the export must precede the run."""
    for path in _scripts_running_the_package():
        lines = path.read_text().splitlines()
        cd_idx = [i for i, line in enumerate(lines) if _CD_TO_ROOT.match(line)]
        assert len(cd_idx) == 1, f"{path.name}: expected one cd-to-repo-root line, got {cd_idx}"
        export_idx = [i for i, line in enumerate(lines) if line == PYTHONPATH_LINE]
        assert len(export_idx) == 1, f"{path.name}: expected one PYTHONPATH export"
        run_idx = [i for i, line in enumerate(lines) if "python -m tps_eval" in line]
        assert cd_idx[0] < export_idx[0] < min(run_idx), (
            f"{path.name}: order must be cd -> export PYTHONPATH -> python -m tps_eval, "
            f"got {cd_idx[0]}, {export_idx[0]}, {min(run_idx)}"
        )
