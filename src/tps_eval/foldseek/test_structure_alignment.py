from __future__ import annotations

"""Unit tests for the structure-alignment tool (foldseek structure_alignment).

Run: python test_structure_alignment.py   (no pytest / foldseek binary required).

foldseek is never executed: subprocess.Popen is monkeypatched with a fake that
writes a tiny synthetic foldseek TSV (the 16-column easy-search format) to the
tool's output path and creates the tmp dir it later removes. That drives main()
end-to-end so the REAL TSV parser, per-query best-hit (idxmax over alntmscore /
qtmscore / ttmscore / lddt) reduction, and the --exclude_self leave-one-out filter
are all exercised on synthetic data. Also unit-tests the pure _structure_stem.
"""

import io
import os
import sys
import tempfile
from types import SimpleNamespace

import pandas as pd


import tps_eval.foldseek.structure_alignment as sa  # noqa: E402
from tps_eval.foldseek.structure_alignment import _structure_stem, main  # noqa: E402

_COLS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart",
    "qend", "tstart", "tend", "evalue", "bits", "alntmscore", "qtmscore",
    "ttmscore", "lddt",
]


class _FakePopen:
    """Writes `tsv_content` to the out_tsv (cmd[4]) and creates tmp dir (cmd[5])."""

    tsv_content = ""

    def __init__(self, cmd, **kwargs):
        out_tsv = cmd[4]
        tmp_dir = cmd[5]
        os.makedirs(tmp_dir, exist_ok=True)
        with open(out_tsv, "w") as fh:
            fh.write(_FakePopen.tsv_content)
        self.stdout = io.StringIO("")

    def wait(self):
        return 0


def _install_fake(tsv_content: str):
    _FakePopen.tsv_content = tsv_content
    sa.subprocess.Popen = _FakePopen  # type: ignore[attr-defined]


def _row(query, target, alntm, qtm, ttm, lddt):
    # Fill the non-tested numeric columns with placeholders.
    return "\t".join(str(v) for v in [
        query, target, 50.0, 100, 10, 1, 1, 100, 1, 100, 1e-20, 200.0,
        alntm, qtm, ttm, lddt,
    ])


def test_structure_stem():
    assert _structure_stem("d1.pdb") == "d1"
    assert _structure_stem("/x/d2.cif") == "d2"
    assert _structure_stem("d3.pdb.gz") == "d3"
    assert _structure_stem("d4.ent") == "d4"
    assert _structure_stem("d5") == "d5"
    print("ok _structure_stem")


def test_best_hit_reduction():
    """Per query, the reported max_* target is the argmax over each score column."""
    tmp = tempfile.mkdtemp(prefix="sa_best_")
    out = os.path.join(tmp, "out")
    tsv = "\n".join([
        # query d1: ref_b wins alntmscore (0.9) and lddt (0.95); ref_a wins qtmscore.
        _row("d1.pdb", "ref_a.pdb", 0.70, 0.99, 0.40, 0.50),
        _row("d1.pdb", "ref_b.pdb", 0.90, 0.60, 0.80, 0.95),
    ]) + "\n"
    _install_fake(tsv)

    args = SimpleNamespace(
        random_run_id=False, output_root=out, structures_root="/fake/q",
        known_structures_root="/fake/k", store_intermediate_results=False,
        exclude_self=False,
    )
    main(args)

    scores = pd.read_csv(os.path.join(out, "structure_alignment_scores.csv")).set_index("query")
    r = scores.loc["d1.pdb"]
    assert r["max_alntmscore"] == 0.90 and r["max_alntmscore_target"] == "ref_b.pdb", r.to_dict()
    assert r["max_qtmscore"] == 0.99 and r["max_qtmscore_target"] == "ref_a.pdb", r.to_dict()
    assert r["max_lddt"] == 0.95 and r["max_lddt_target"] == "ref_b.pdb", r.to_dict()
    print("ok best-hit reduction (per-column argmax)")


def test_exclude_self_leave_one_out():
    """--exclude_self drops target-stem == query-stem before the best-hit reduction,
    so a set searched against itself returns the nearest OTHER neighbour."""
    tmp = tempfile.mkdtemp(prefix="sa_self_")
    out = os.path.join(tmp, "out")
    tsv = "\n".join([
        # self-hit (TM 1.0) must be dropped; the real neighbour ref9 (0.72) wins.
        _row("d1.pdb", "d1.pdb", 1.0, 1.0, 1.0, 1.0),
        _row("d1.pdb", "ref9.pdb", 0.72, 0.70, 0.71, 0.68),
    ]) + "\n"
    _install_fake(tsv)

    args = SimpleNamespace(
        random_run_id=False, output_root=out, structures_root="/fake/q",
        known_structures_root="/fake/k", store_intermediate_results=False,
        exclude_self=True,
    )
    main(args)

    scores = pd.read_csv(os.path.join(out, "structure_alignment_scores.csv")).set_index("query")
    r = scores.loc["d1.pdb"]
    assert r["max_alntmscore"] == 0.72, r.to_dict()          # NOT the 1.0 self-hit
    assert r["max_alntmscore_target"] == "ref9.pdb", r.to_dict()
    print("ok exclude_self leave-one-out")


def test_all_nan_score_column_yields_nan_not_crash():
    """A query whose entire lddt column is non-numeric (-> all NaN after coercion) yields a
    NaN idxmax; the tool must report NaN max_lddt for it, not crash on df.loc[NaN] (regression)."""
    tmp = tempfile.mkdtemp(prefix="sa_nan_")
    out = os.path.join(tmp, "out")
    tsv = "\n".join([
        _row("d1.pdb", "ref_a.pdb", 0.70, 0.60, 0.50, "NA"),
        _row("d1.pdb", "ref_b.pdb", 0.90, 0.80, 0.70, "NA"),
    ]) + "\n"
    _install_fake(tsv)

    args = SimpleNamespace(
        random_run_id=False, output_root=out, structures_root="/fake/q",
        known_structures_root="/fake/k", store_intermediate_results=False,
        exclude_self=False,
    )
    main(args)  # must not raise

    scores = pd.read_csv(os.path.join(out, "structure_alignment_scores.csv")).set_index("query")
    r = scores.loc["d1.pdb"]
    assert pd.isna(r["max_lddt"]), r.to_dict()        # all-NaN column -> NaN, no KeyError
    assert r["max_alntmscore"] == 0.90, r.to_dict()   # numeric columns still reduce
    print("ok all-NaN score column -> NaN (no crash)")


def main_all():
    test_structure_stem()
    test_best_hit_reduction()
    test_exclude_self_leave_one_out()
    test_all_nan_score_column_yields_nan_not_crash()
    print("\nAll 4 tests passed.")


if __name__ == "__main__":
    main_all()
