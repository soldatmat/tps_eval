"""The src-layout anchors in ``repo_paths`` must point at things that exist.

Regression: after the src-layout migration several modules computed the repo root
as ``Path(__file__).parent.parent.parent`` — one level short — and silently resolved
to ``<repo>/src``. Nothing raised: ProteinMPNN was invoked as
``<repo>/src/vendor/ProteinMPNN/protein_mpnn_run.py``, the subprocess failed per
structure, the failure was caught per-design and the tool wrote an ALL-NaN column
that looked like a legitimate "no signal" result (two production runs shipped that
way). These assertions are the cheap alarm that would have caught it.
"""
from __future__ import annotations

from pathlib import Path

from tps_eval import repo_paths


def test_anchors_are_nested_as_documented():
    assert repo_paths.PACKAGE_DIR.name == "tps_eval"
    assert repo_paths.SRC_DIR.name == "src"
    assert repo_paths.PACKAGE_DIR.parent == repo_paths.SRC_DIR
    assert repo_paths.SRC_DIR.parent == repo_paths.REPO_ROOT
    # The repo root is a checkout, not the source root.
    assert (repo_paths.REPO_ROOT / "pyproject.toml").is_file()
    assert not (repo_paths.REPO_ROOT / "tps_eval").exists()


def test_repo_level_dirs_exist():
    assert repo_paths.SCRIPTS_DIR.is_dir()
    assert (repo_paths.SCRIPTS_DIR / "submit_job.sh").is_file()
    # vendor/ holds git submodules; the dir exists even before `submodule update`.
    assert repo_paths.VENDOR_DIR.is_dir()


def test_package_data_dirs_exist():
    assert repo_paths.REFERENCE_STATS_DIR.is_dir()
    assert list(repo_paths.REFERENCE_STATS_DIR.glob("marts_db_*_metric_stats.json"))


def test_vendor_shellouts_resolve_under_the_repo_not_src():
    """Every module that shells out to a vendored tool must escape src/."""
    from tps_eval.sequence_metrics.catapro import CATAPRO_DIR, CATAPRO_PREDICT
    from tps_eval.structure_metrics.proteinmpnn_score import PROTEINMPNN_DIR
    from tps_eval.structure_metrics.self_consistency import (
        PROTEINMPNN_DIR as SC_PROTEINMPNN_DIR,
    )

    for path in (PROTEINMPNN_DIR, SC_PROTEINMPNN_DIR, CATAPRO_DIR, CATAPRO_PREDICT):
        assert repo_paths.SRC_DIR not in Path(path).parents, (
            f"{path} resolves inside src/ — the repo root was computed one level short"
        )
    assert PROTEINMPNN_DIR == repo_paths.VENDOR_DIR / "ProteinMPNN"
    assert SC_PROTEINMPNN_DIR == repo_paths.VENDOR_DIR / "ProteinMPNN"
    assert CATAPRO_DIR == repo_paths.VENDOR_DIR / "CataPro"


def test_dashboard_finds_its_default_reference_bands():
    from tps_eval.dashboard import build_dashboard

    bands = build_dashboard.default_source_paths()
    assert bands, "dashboard would render with NO MARTS-DB reference bands"
    assert all(Path(p).is_file() for p in bands)
