from __future__ import annotations

"""Self-contained tests for max_sequence_identity.py.

Run from this directory (flat-module import resolves like the runner does; the
module inserts src/ onto sys.path for its `from data.sequences` import):
    cd src/sequence_metrics && python test_max_sequence_identity.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_max_sequence_identity.py -q

No network / no models — uses Biopython's bundled BLOSUM62 (loaded at import) and
tiny synthetic sequences. Locks in: global-alignment identity/similarity on known
pairs (identical -> 1.0/1.0; all A vs all S -> identity 0 but similarity 1 since
BLOSUM62 A/S is positive), the two-set max reduction + hit index, self mode
excluding the same index, empty-generated robustness, the top-k percent-score
contract, ID keying + CSV output filename, and the save-path helpers.

Note: the reduction runs a real (small) ProcessPoolExecutor fork pool, so it
exercises the actual parallel code path.
"""

import os
import tempfile

from max_sequence_identity import (
    _get_save_path,
    _get_topk_save_path,
    _pair_metrics,
    evaluate_max_sequence_identity,
    get_max_sequence_identity,
    get_max_sequence_identity_two_sets,
    get_topk_sequence_identity,
)

import pandas as pd


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b}"


def test_pair_metrics_identical():
    ident, sim = _pair_metrics("ACDEFGHIK", "ACDEFGHIK")
    _approx(ident, 1.0)
    _approx(sim, 1.0)


def test_pair_metrics_positive_but_not_identical():
    # BLOSUM62 A/S == +1 -> similarity counts, identity does not.
    ident, sim = _pair_metrics("AAAA", "SSSS")
    _approx(ident, 0.0)
    _approx(sim, 1.0)


def test_two_sets_picks_best_and_hit_index():
    train = ["MMMMMMMMM", "ACDEFGHIK"]  # index 1 is the identical match
    generated = ["ACDEFGHIK"]
    ident, sim, ident_idx, sim_idx = get_max_sequence_identity_two_sets(train, generated)
    _approx(ident[0], 1.0)
    assert ident_idx[0] == 1  # best identity hit is train[1]


def test_self_mode_excludes_same_index():
    # In self mode a sequence must not match itself; here two entries are identical
    # so each still finds a perfect *other* match at the twin index.
    seqs = ["ACDEFGHIK", "ACDEFGHIK", "MMMMMMMMM"]
    ident, _sim, ident_idx, _si = get_max_sequence_identity(seqs, self_comparison=True)
    _approx(ident[0], 1.0)
    assert ident_idx[0] == 1  # matched its twin, not itself (index 0)
    assert ident_idx[1] == 0


def test_empty_generated_returns_empty_lists():
    out = get_max_sequence_identity_two_sets(["ACDE"], [])
    assert out == ([], [], [], [])


def test_topk_percent_scores_and_order():
    train = ["ACDEFGHIK", "ACDEFGHIL", "MMMMMMMMM"]
    generated = ["ACDEFGHIK"]
    topk = get_topk_sequence_identity(train, generated, top_k=2)
    ranked = topk[0]
    assert len(ranked) == 2
    # Scores are PERCENT; the exact identical hit (train[0]) is first at 100%.
    (idx0, score0), (idx1, score1) = ranked
    assert idx0 == 0
    _approx(score0, 100.0, tol=1e-6)
    assert score0 >= score1


def test_evaluate_writes_csv_keyed_by_id():
    train = ["MMMMMMMMM", "ACDEFGHIK"]
    generated = ["ACDEFGHIK", "MMMMMMMMM"]
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "out.csv")
        evaluate_max_sequence_identity(
            train, generated,
            generated_identifiers=["g1", "g2"],
            train_identifiers=["t0", "t1"],
            save_path=sp,
        )
        assert os.path.isfile(sp)
        df = pd.read_csv(sp).set_index("ID")
        assert set(df.columns) >= {
            "sequence_identity", "sequence_identity_hit",
            "sequence_similarity", "sequence_similarity_hit",
        }
        _approx(float(df.loc["g1", "sequence_identity"]), 1.0)
        assert df.loc["g1", "sequence_identity_hit"] == "t1"  # train[1] == ACDEFGHIK
        assert df.loc["g2", "sequence_identity_hit"] == "t0"  # train[0] == MMMMMMMMM


def test_topk_csv_written():
    train = ["ACDEFGHIK", "MMMMMMMMM"]
    generated = ["ACDEFGHIK"]
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "out.csv")
        tp = os.path.join(d, "out_topk.csv")
        evaluate_max_sequence_identity(
            train, generated,
            generated_identifiers=["g1"],
            train_identifiers=["t0", "t1"],
            save_path=sp,
            top_k=2,
            topk_save_path=tp,
        )
        assert os.path.isfile(tp)
        topk = pd.read_csv(tp)
        assert list(topk.columns) == ["query_id", "rank", "neighbour_id", "score"]
        assert list(topk["query_id"]) == ["g1", "g1"]
        # Best hit is the identical train[0] at rank 1 with 100% score.
        assert topk.iloc[0]["neighbour_id"] == "t0"
        _approx(float(topk.iloc[0]["score"]), 100.0, tol=1e-6)


def test_save_path_naming():
    assert _get_save_path("designs.fasta") == "designs_max_sequence_identity.csv"
    assert _get_save_path("designs.fasta", save_suffix="self") == (
        "designs_max_sequence_identity_self.csv"
    )
    assert _get_topk_save_path("designs.fasta") == (
        "designs_max_sequence_identity_topk.csv"
    )


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
