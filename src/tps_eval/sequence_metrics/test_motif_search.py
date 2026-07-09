from __future__ import annotations

"""Self-contained tests for motif_search.py.

Run from this directory (flat-module import resolves like the runner does; the
module itself inserts src/ onto sys.path for its `from data.sequences` import):
    cd src/sequence_metrics && python test_motif_search.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_motif_search.py -q

Locks in: per-motif boolean columns keyed by the regex pattern string, ID keying
from the FASTA, both string-motif and pre-compiled-pattern inputs, the
`<input>_motifs.csv` save path + that the file is written, and empty-FASTA
robustness. Uses only tiny temp FASTA files (no network, no models).
"""

import os
import re
import tempfile

from tps_eval.sequence_metrics.motif_search import convert_motifs, find_motifs, get_save_path, motif_search


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for ident, seq in records:
            fh.write(f">{ident}\n{seq}\n")


def test_get_save_path():
    assert get_save_path("designs.fasta") == "designs_motifs.csv"
    assert get_save_path("/a/b/designs.fa") == "/a/b/designs_motifs.csv"


def test_convert_motifs_compiles():
    patterns = convert_motifs(["DD..D", "M+"])
    assert all(isinstance(p, re.Pattern) for p in patterns)
    assert patterns[0].pattern == "DD..D"


def test_string_motifs_boolean_columns():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "seqs.fasta")
        _write_fasta(p, [("s1", "AADDLLDAA"), ("s2", "GGGGGGGG"), ("s3", "MMMK")])
        df = motif_search(p, ["DD..D", "M+"], save=True)
        assert list(df["ID"]) == ["s1", "s2", "s3"]
        # Columns are named by the regex pattern string.
        assert "DD..D" in df.columns and "M+" in df.columns
        by_id = df.set_index("ID")
        assert bool(by_id.loc["s1", "DD..D"]) is True
        assert bool(by_id.loc["s2", "DD..D"]) is False
        assert bool(by_id.loc["s1", "M+"]) is False
        assert bool(by_id.loc["s3", "M+"]) is True
        # CSV written next to the input.
        assert os.path.isfile(os.path.join(d, "seqs_motifs.csv"))


def test_precompiled_pattern_input():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "seqs.fasta")
        _write_fasta(p, [("s1", "AADDLLDAA")])
        df = motif_search(p, convert_motifs(["DD..D"]), save=False)
        assert bool(df.set_index("ID").loc["s1", "DD..D"]) is True
        # save=False -> no CSV written.
        assert not os.path.isfile(os.path.join(d, "seqs_motifs.csv"))


def test_find_motifs_in_memory():
    import pandas as pd

    frame = pd.DataFrame({"ID": ["a", "b"], "sequence": ["DDLLD", "AAAA"]})
    find_motifs(frame, convert_motifs(["DD..D"]))
    assert list(frame["DD..D"]) == [True, False]


def test_empty_fasta_no_crash():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "empty.fasta")
        open(p, "w").close()
        df = motif_search(p, ["DD..D"], save=True)
        assert list(df["ID"]) == []
        assert "DD..D" in df.columns
        assert os.path.isfile(os.path.join(d, "empty_motifs.csv"))


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
