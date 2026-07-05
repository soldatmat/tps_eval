from __future__ import annotations

"""Unit tests for the broad STRUCTURE homology search tool (foldseek_swissprot_search).

Run: python test_foldseek_swissprot_search.py   (no pytest / foldseek binary required).

foldseek is never executed: subprocess.Popen is monkeypatched with a fake that
records the built command and writes a tiny synthetic foldseek TSV to the tool's
output path. Exercises the REAL command construction, AFDB-accession extraction,
query->ID mapping, best-hit-by-alntmscore selection, TPS classification, ID keying,
empty-result NaN contract, and CSV filename — all on synthetic data.
"""

import io
import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # src/

from homology_search import foldseek_swissprot_search as fss  # noqa: E402
from homology_search.foldseek_swissprot_search import (  # noqa: E402
    _accession_from_target,
    _default_save_path,
    _query_to_id,
    _summarize_query,
    foldseek_swissprot_search,
)


class _FakePopen:
    captured_cmd = None
    tsv_content = ""

    def __init__(self, cmd, **kwargs):
        _FakePopen.captured_cmd = list(cmd)
        # foldseek easy-search <query> <db> <out_tsv> <tmp> ...  -> out is index 4.
        out_path = cmd[4]
        with open(out_path, "w") as fh:
            fh.write(_FakePopen.tsv_content)
        self.stdout = io.StringIO("")

    def wait(self):
        return 0


def _install_fake(tsv_content: str):
    _FakePopen.tsv_content = tsv_content
    fss.subprocess.Popen = _FakePopen  # type: ignore[attr-defined]


def _write_pdb(path: str) -> None:
    """A minimal one-CA .pdb (never parsed by the tool — _collect_structures only
    globs; foldseek is faked — so contents just need to be a valid-ish file)."""
    with open(path, "w") as fh:
        fh.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")


def _write_accessions(lines) -> str:
    fd, path = tempfile.mkstemp(prefix="fss_acc_", suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------- #
# Pure helpers.
# --------------------------------------------------------------------------- #
def test_accession_from_target():
    assert _accession_from_target("AF-P12345-F1-model_v4") == "P12345"
    assert _accession_from_target("AF-Q9ABC1-F2-model_v3.pdb") == "Q9ABC1"
    assert _accession_from_target("plain.cif") == "plain"       # non-AFDB -> stem
    assert _accession_from_target("noext") == "noext"
    print("ok _accession_from_target")


def test_query_to_id():
    assert _query_to_id("d1.pdb") == "d1"
    assert _query_to_id("/tmp/q/d2.cif") == "d2"
    assert _query_to_id("job7_model.cif") == "job7"            # AF3 model stem
    assert _query_to_id("d3") == "d3"
    print("ok _query_to_id")


def test_default_save_path():
    assert _default_save_path("/x/gen_structs") == "/x/gen_structs_foldseek_swissprot_search.csv"
    assert _default_save_path("/x/gen_structs/") == "/x/gen_structs_foldseek_swissprot_search.csv"
    print("ok _default_save_path")


def test_summarize_query_tmscore_ranking():
    """Best hit = MAX alntmscore; best_nontps_tmscore + TPS count."""
    tps_set = frozenset({"P0C2A9"})
    group = pd.DataFrame(
        {
            "target": ["AF-Q00001-F1-model_v4", "AF-P0C2A9-F1-model_v4", "AF-Q00002-F1-model_v4"],
            "alntmscore": [0.60, 0.85, 0.40],
            "bits": [200, 300, 100],
        }
    )
    row = _summarize_query(group, tps_set, top_n=25)
    assert row["foldseek_sprot_top_hit"] == "P0C2A9", row      # highest TM wins
    assert row["foldseek_sprot_top_tmscore"] == 0.85, row
    assert row["foldseek_sprot_top_is_tps"] is True, row
    assert row["foldseek_sprot_best_nontps_tmscore"] == 0.60, row
    assert row["foldseek_sprot_n_tps_in_topN"] == 1, row
    print("ok _summarize_query TM ranking + classification")


# --------------------------------------------------------------------------- #
# Full run through the (faked) foldseek call.
# --------------------------------------------------------------------------- #
def test_full_run_command_and_parsing():
    tmp = tempfile.mkdtemp(prefix="fss_full_")
    structs = os.path.join(tmp, "gen_structs")
    os.makedirs(structs)
    _write_pdb(os.path.join(structs, "d1.pdb"))
    _write_pdb(os.path.join(structs, "d2.pdb"))  # d2 will have no hits
    acc = _write_accessions(["P0C2A9"])

    tsv = (
        "d1.pdb\tAF-P0C2A9-F1-model_v4\t0.8\t200\t1e-50\t300\t0.85\n"
        "d1.pdb\tAF-Q00001-F1-model_v4\t0.5\t180\t1e-30\t200\t0.60\n"
    )
    _install_fake(tsv)

    df = foldseek_swissprot_search(structs, "/fake/afdb", acc, top_n=25, max_seqs=300)

    # --- command construction ---
    cmd = _FakePopen.captured_cmd
    assert cmd[0:2] == ["foldseek", "easy-search"], cmd
    assert cmd[cmd.index("--max-seqs") + 1] == "300", cmd
    assert cmd[cmd.index("--format-output") + 1] == fss.FOLDSEEK_OUTFMT, cmd

    # --- parsing / classification ---
    by_id = df.set_index("ID")
    assert set(df["ID"]) == {"d1", "d2"}
    assert by_id.loc["d1", "foldseek_sprot_top_hit"] == "P0C2A9"
    assert by_id.loc["d1", "foldseek_sprot_top_tmscore"] == 0.85
    assert bool(by_id.loc["d1", "foldseek_sprot_top_is_tps"]) is True
    assert by_id.loc["d1", "foldseek_sprot_best_nontps_tmscore"] == 0.60
    assert int(by_id.loc["d1", "foldseek_sprot_n_tps_in_topN"]) == 1
    # d2: no hits.
    assert by_id.loc["d2", "foldseek_sprot_top_hit"] == "" or pd.isna(by_id.loc["d2", "foldseek_sprot_top_hit"])
    assert int(by_id.loc["d2", "foldseek_sprot_n_tps_in_topN"]) == 0
    assert pd.isna(by_id.loc["d2", "foldseek_sprot_top_tmscore"])

    assert os.path.exists(_default_save_path(structs))
    print("ok full run: command + parsing + NaN-on-no-hit + CSV")


def test_full_run_empty_results():
    tmp = tempfile.mkdtemp(prefix="fss_empty_")
    structs = os.path.join(tmp, "gen_structs")
    os.makedirs(structs)
    _write_pdb(os.path.join(structs, "d1.pdb"))
    acc = _write_accessions(["P0C2A9"])
    _install_fake("")

    df = foldseek_swissprot_search(structs, "/fake/afdb", acc)
    assert len(df) == 1
    r = df.set_index("ID").loc["d1"]
    assert int(r["foldseek_sprot_n_tps_in_topN"]) == 0
    assert pd.isna(r["foldseek_sprot_top_tmscore"])
    print("ok full run: empty results -> NaN not crash")


def main():
    test_accession_from_target()
    test_query_to_id()
    test_default_save_path()
    test_summarize_query_tmscore_ranking()
    test_full_run_command_and_parsing()
    test_full_run_empty_results()
    print("\nAll 6 tests passed.")


if __name__ == "__main__":
    main()
