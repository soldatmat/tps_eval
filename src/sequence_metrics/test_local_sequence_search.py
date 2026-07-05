from __future__ import annotations

"""Self-contained tests for local_sequence_search.py.

Run from this directory (flat-module import resolves like the runner does; the
module inserts src/ onto sys.path for its `from data.sequences` import):
    cd src/sequence_metrics && python test_local_sequence_search.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_local_sequence_search.py -q

The mmseqs2/diamond BACKENDS are external binaries and are NOT invoked here (they
are absent on this machine) — instead these tests exercise the pure-python
reduction/parse layer that the backends feed into, which is where the real
result-shaping logic lives: best-hit-per-query selection (max bitscore, self
exclusion, deterministic tie-break), the tidy top-k neighbours table (dedup +
score-descending ranking), first-token ID normalization, the unknown-backend
guard, empty-hits robustness, and the `<input>_local_sequence_search*.csv` save
paths. All inputs are synthetic in-memory hit DataFrames.
"""

import numpy as np
import pandas as pd

from local_sequence_search import (
    BACKENDS,
    TOPK_COLUMNS,
    _HIT_COLUMNS,
    _best_hits,
    _default_save_path,
    _default_topk_save_path,
    _first_token,
    _topk_neighbours,
    local_sequence_search,
)


def _hits(rows):
    """Build a hits frame with the internal normalized columns (+ _bits)."""
    cols = ["qseqid", "sseqid", "identity", "similarity", "coverage", "score", "_bits"]
    return pd.DataFrame(rows, columns=cols)


def test_first_token():
    assert _first_token("seq1 some description") == "seq1"
    assert _first_token("seq2\ttab") == "seq2"
    assert _first_token("plain") == "plain"


def test_unknown_backend_raises():
    try:
        local_sequence_search("x.fasta", backend="nope")
    except ValueError as e:
        assert "nope" in str(e)
    else:
        raise AssertionError("expected ValueError for unknown backend")
    assert set(BACKENDS) == {"mmseqs2", "diamond"}


def test_best_hits_picks_max_bitscore():
    hits = _hits([
        ["q1", "a", 40.0, 60.0, 90.0, 40.0, 100.0],
        ["q1", "b", 80.0, 90.0, 95.0, 80.0, 200.0],  # higher bits -> best
        ["q2", "c", 55.0, 70.0, 88.0, 55.0, 150.0],
    ])
    best = _best_hits(hits, self_mode=False)
    by_q = best.set_index("qseqid")
    assert by_q.loc["q1", "sseqid"] == "b"
    assert by_q.loc["q1", "identity"] == 80.0
    assert by_q.loc["q2", "sseqid"] == "c"


def test_best_hits_self_mode_excludes_self():
    hits = _hits([
        ["q1", "q1", 100.0, 100.0, 100.0, 100.0, 999.0],  # self -> excluded
        ["q1", "b", 70.0, 80.0, 90.0, 70.0, 120.0],
    ])
    best = _best_hits(hits, self_mode=True)
    assert list(best["sseqid"]) == ["b"]


def test_best_hits_empty():
    empty = pd.DataFrame(columns=_HIT_COLUMNS + ["_bits"])
    assert _best_hits(empty, self_mode=False).empty


def test_topk_neighbours_ranking_and_dedup():
    hits = _hits([
        ["q1", "a", 90.0, 0.0, 0.0, 90.0, 300.0],
        ["q1", "a", 90.0, 0.0, 0.0, 90.0, 250.0],  # duplicate (q1,a) -> keep best bits
        ["q1", "b", 70.0, 0.0, 0.0, 70.0, 200.0],
        ["q1", "c", 50.0, 0.0, 0.0, 50.0, 100.0],
    ])
    topk = _topk_neighbours(hits, ["q1"], top_k=2, self_mode=False)
    assert list(topk.columns) == TOPK_COLUMNS
    # Score (identity %) descending: a (90) rank1, b (70) rank2; c dropped by top_k.
    assert list(topk["neighbour_id"]) == ["a", "b"]
    assert list(topk["rank"]) == [1, 2]
    np.testing.assert_allclose(topk["score"], [90.0, 70.0])


def test_topk_neighbours_self_mode():
    hits = _hits([
        ["q1", "q1", 100.0, 0.0, 0.0, 100.0, 999.0],  # self -> excluded
        ["q1", "b", 60.0, 0.0, 0.0, 60.0, 120.0],
    ])
    topk = _topk_neighbours(hits, ["q1"], top_k=5, self_mode=True)
    assert list(topk["neighbour_id"]) == ["b"]


def test_topk_neighbours_empty():
    empty = pd.DataFrame(columns=_HIT_COLUMNS + ["_bits"])
    out = _topk_neighbours(empty, ["q1", "q2"], top_k=3, self_mode=False)
    assert list(out.columns) == TOPK_COLUMNS
    assert len(out) == 0


def test_save_path_naming():
    assert _default_save_path("designs.fasta") == "designs_local_sequence_search.csv"
    assert (
        _default_topk_save_path("designs.fasta")
        == "designs_local_sequence_search_topk.csv"
    )


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
