from __future__ import annotations

"""Unit tests for the broad SEQUENCE homology search tool (swissprot_search).

Run: python test_swissprot_search.py   (no pytest / DIAMOND binary required).

DIAMOND is never executed: subprocess.Popen is monkeypatched with a fake that
records the built command and writes a tiny synthetic blast-tab (outfmt 6) file to
the tool's --out path. That exercises the REAL command construction, m8/TSV parser,
best-hit-by-bitscore selection, TPS/non-TPS classification, ID keying, empty-result
NaN contract, and CSV filename — all on synthetic data.
"""

import io
import os
import sys
import tempfile

import pandas as pd


from tps_eval.homology_search import swissprot_search as sws  # noqa: E402
from tps_eval.homology_search.swissprot_search import (  # noqa: E402
    _accession_from_sseqid,
    _default_save_path,
    _summarize_query,
    swissprot_search,
)


# --------------------------------------------------------------------------- #
# Fake DIAMOND: captures the command + writes synthetic hits to the --out path.
# --------------------------------------------------------------------------- #
class _FakePopen:
    captured_cmd = None
    tsv_content = ""  # set per-test

    def __init__(self, cmd, **kwargs):
        _FakePopen.captured_cmd = list(cmd)
        out_path = cmd[cmd.index("--out") + 1]
        with open(out_path, "w") as fh:
            fh.write(_FakePopen.tsv_content)
        self.stdout = io.StringIO("")

    def wait(self):
        return 0


def _install_fake(tsv_content: str):
    _FakePopen.tsv_content = tsv_content
    sws.subprocess.Popen = _FakePopen  # type: ignore[attr-defined]


def _write_fasta(path: str, records) -> None:
    with open(path, "w") as fh:
        for ident, seq in records:
            fh.write(f">{ident}\n{seq}\n")


def _write_accessions(lines) -> str:
    fd, path = tempfile.mkstemp(prefix="sws_acc_", suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------- #
# Pure helpers.
# --------------------------------------------------------------------------- #
def test_accession_from_sseqid():
    assert _accession_from_sseqid("sp|P12345|NAME_SP") == "P12345"
    assert _accession_from_sseqid("tr|Q9ABC1|OTHER") == "Q9ABC1"
    assert _accession_from_sseqid("P12345") == "P12345"          # plain-accession DB
    assert _accession_from_sseqid("weird|token") == "weird|token"  # <3 fields -> verbatim
    print("ok _accession_from_sseqid")


def test_default_save_path():
    assert _default_save_path("/x/designs.fasta") == "/x/designs_swissprot_search.csv"
    print("ok _default_save_path")


def test_summarize_query_bitscore_and_classification():
    """Best hit = MAX bitscore (not max pident); best_nontps_pident + TPS count."""
    tps_set = frozenset({"P0C2A9"})
    # A non-TPS hit has the highest pident (95) but a LOWER bitscore than the TPS hit.
    group = pd.DataFrame(
        {
            "qseqid": ["d1", "d1", "d1"],
            "sseqid": ["sp|Q00001|NONTPS", "sp|P0C2A9|TPS", "sp|Q00002|NONTPS2"],
            "pident": [95.0, 80.0, 40.0],
            "bitscore": [200.0, 300.0, 100.0],
            "evalue": [1e-40, 1e-90, 1e-10],
        }
    )
    row = _summarize_query(group, tps_set, top_n=25)
    assert row["swissprot_top_hit"] == "P0C2A9", row       # highest bitscore wins
    assert row["swissprot_top_pident"] == 80.0, row
    assert row["swissprot_top_bitscore"] == 300.0, row
    assert row["swissprot_top_is_tps"] is True, row
    assert row["swissprot_best_nontps_pident"] == 95.0, row  # best pident among non-TPS
    assert row["swissprot_n_tps_in_topN"] == 1, row
    print("ok _summarize_query best-hit + classification")


def test_summarize_query_all_nontps():
    tps_set = frozenset({"P0C2A9"})
    group = pd.DataFrame(
        {
            "qseqid": ["d2", "d2"],
            "sseqid": ["sp|Q1|A", "sp|Q2|B"],
            "pident": [50.0, 70.0],
            "bitscore": [120.0, 90.0],
            "evalue": [1e-20, 1e-15],
        }
    )
    row = _summarize_query(group, tps_set, top_n=25)
    assert row["swissprot_top_is_tps"] is False, row
    assert row["swissprot_n_tps_in_topN"] == 0, row
    assert row["swissprot_best_nontps_pident"] == 70.0, row
    print("ok _summarize_query all-non-TPS")


# --------------------------------------------------------------------------- #
# Full run through the (faked) DIAMOND call.
# --------------------------------------------------------------------------- #
def test_full_run_command_and_parsing():
    tmp = tempfile.mkdtemp(prefix="sws_full_")
    fasta = os.path.join(tmp, "designs.fasta")
    # d1 has hits; d2 has NO hits (must still get a NaN/empty row).
    _write_fasta(fasta, [("d1", "MTTYVKLANDE"), ("d2", "MSSSAAAAAAA")])
    acc = _write_accessions(["P0C2A9"])

    tsv = (
        "d1\tsp|P0C2A9|TPS_SP\t85.0\t300.0\t1e-90\n"
        "d1\tsp|Q00001|OTHER_SP\t60.0\t200.0\t1e-40\n"
    )
    _install_fake(tsv)

    df = swissprot_search(
        fasta, "/fake/db", acc, top_n=25, threads=2, sensitivity="very-sensitive"
    )

    # --- command construction ---
    cmd = _FakePopen.captured_cmd
    assert cmd[0:2] == ["diamond", "blastp"], cmd
    assert "--very-sensitive" in cmd, cmd
    assert cmd[cmd.index("--max-target-seqs") + 1] == "25", cmd
    assert cmd[cmd.index("--threads") + 1] == "2", cmd
    # outfmt 6 followed by the exact column spec the parser relies on.
    fmt_i = cmd.index("--outfmt")
    assert cmd[fmt_i + 1] == "6", cmd
    assert cmd[fmt_i + 2 : fmt_i + 2 + len(sws.DIAMOND_OUTFMT)] == sws.DIAMOND_OUTFMT, cmd

    # --- parsing / classification ---
    by_id = df.set_index("ID")
    assert set(df["ID"]) == {"d1", "d2"}
    assert by_id.loc["d1", "swissprot_top_hit"] == "P0C2A9"
    assert by_id.loc["d1", "swissprot_top_bitscore"] == 300.0
    assert bool(by_id.loc["d1", "swissprot_top_is_tps"]) is True
    assert by_id.loc["d1", "swissprot_best_nontps_pident"] == 60.0
    assert int(by_id.loc["d1", "swissprot_n_tps_in_topN"]) == 1
    # d2: no hits -> empty top hit, 0 TPS in topN, NA is_tps.
    assert by_id.loc["d2", "swissprot_top_hit"] == "" or pd.isna(by_id.loc["d2", "swissprot_top_hit"])
    assert int(by_id.loc["d2", "swissprot_n_tps_in_topN"]) == 0
    assert pd.isna(by_id.loc["d2", "swissprot_top_pident"])

    # --- CSV written to the default sibling path ---
    assert os.path.exists(_default_save_path(fasta))
    print("ok full run: command + parsing + NaN-on-no-hit + CSV")


def test_full_run_empty_results():
    """An empty DIAMOND output must yield all-NaN rows, not a crash."""
    tmp = tempfile.mkdtemp(prefix="sws_empty_")
    fasta = os.path.join(tmp, "designs.fasta")
    _write_fasta(fasta, [("d1", "MTTYVKLANDE")])
    acc = _write_accessions(["P0C2A9"])
    _install_fake("")  # zero-byte output

    df = swissprot_search(fasta, "/fake/db", acc)
    assert len(df) == 1
    r = df.set_index("ID").loc["d1"]
    assert int(r["swissprot_n_tps_in_topN"]) == 0
    assert pd.isna(r["swissprot_top_bitscore"])
    print("ok full run: empty results -> NaN not crash")


def main():
    test_accession_from_sseqid()
    test_default_save_path()
    test_summarize_query_bitscore_and_classification()
    test_summarize_query_all_nontps()
    test_full_run_command_and_parsing()
    test_full_run_empty_results()
    print("\nAll 6 tests passed.")


if __name__ == "__main__":
    main()
