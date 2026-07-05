from __future__ import annotations

"""Unit tests for the domain-alignment tool (foldseek domain_alignment).

Run: python test_domain_alignment.py   (no pytest / foldseek binary required).

foldseek is never executed: subprocess.Popen is monkeypatched with a fake that
writes a tiny synthetic foldseek TSV to the tool's output path and creates the tmp
dir it later removes, driving main() end-to-end. Exercises the REAL 16-column TSV
parser and the per-query best-hit (idxmax over alntmscore / qtmscore / ttmscore /
lddt) reduction on synthetic data.
"""

import io
import os
import sys
import tempfile
from types import SimpleNamespace

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import domain_alignment as da  # noqa: E402
from domain_alignment import main  # noqa: E402


class _FakePopen:
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
    da.subprocess.Popen = _FakePopen  # type: ignore[attr-defined]


def _row(query, target, alntm, qtm, ttm, lddt):
    return "\t".join(str(v) for v in [
        query, target, 50.0, 100, 10, 1, 1, 100, 1, 100, 1e-20, 200.0,
        alntm, qtm, ttm, lddt,
    ])


def test_best_hit_reduction_and_multiple_queries():
    tmp = tempfile.mkdtemp(prefix="da_best_")
    out = os.path.join(tmp, "out")
    tsv = "\n".join([
        _row("dom1.pdb", "known_alpha.pdb", 0.55, 0.60, 0.50, 0.40),
        _row("dom1.pdb", "known_beta.pdb", 0.88, 0.50, 0.90, 0.92),
        _row("dom2.pdb", "known_alpha.pdb", 0.30, 0.35, 0.33, 0.31),
    ]) + "\n"
    _install_fake(tsv)

    args = SimpleNamespace(
        random_run_id=False, output_root=out,
        detected_domain_structures_root="/fake/detected",
        known_domain_structures_root="/fake/known",
        store_intermediate_results=False,
    )
    main(args)

    scores = pd.read_csv(os.path.join(out, "domain_alignment_scores.csv")).set_index("query")
    assert set(scores.index) == {"dom1.pdb", "dom2.pdb"}, scores.index.tolist()
    r1 = scores.loc["dom1.pdb"]
    assert r1["max_alntmscore"] == 0.88 and r1["max_alntmscore_target"] == "known_beta.pdb", r1.to_dict()
    assert r1["max_qtmscore"] == 0.60 and r1["max_qtmscore_target"] == "known_alpha.pdb", r1.to_dict()
    assert r1["max_lddt"] == 0.92 and r1["max_lddt_target"] == "known_beta.pdb", r1.to_dict()
    # dom2 has a single hit -> that hit is the best across all columns.
    r2 = scores.loc["dom2.pdb"]
    assert r2["max_alntmscore"] == 0.30 and r2["max_alntmscore_target"] == "known_alpha.pdb", r2.to_dict()
    print("ok best-hit reduction + multiple queries")


def main_all():
    test_best_hit_reduction_and_multiple_queries()
    print("\nAll 1 tests passed.")


if __name__ == "__main__":
    main_all()
