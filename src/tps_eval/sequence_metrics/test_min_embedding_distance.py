from __future__ import annotations

"""Self-contained tests for min_embedding_distance.py.

Run from this directory (flat-module import resolves like the runner does; the
module inserts src/ onto sys.path for its `from data.embeddings` import):
    cd src/sequence_metrics && python test_min_embedding_distance.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_min_embedding_distance.py -q

No torch/ESM needed — this tool consumes precomputed embedding CSVs. Tests use
closed-form L2 distances on tiny 2-D synthetic embeddings: pairwise distance
matrix, argmin nearest-neighbour selection + hit ID, self mode (diagonal masked
to +inf so a sequence never matches itself), the tidy top-k neighbours CSV
(ascending distance, self excluded), the `_min_embedding_distance*.csv` save
paths, and a load_embeddings CSV round-trip.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from tps_eval.sequence_metrics.min_embedding_distance import (
    _get_save_path,
    _get_topk_save_path,
    _min_embedding_distance,
    _min_embedding_distance_self,
    get_distances,
    get_min_distances,
    main_train_sequences,
    preprocess_embeddings,
    save_embeddings,
    write_topk_distances,
)


def _emb_df(ids, vectors):
    return pd.DataFrame({"ID": list(ids), "embedding": [list(v) for v in vectors]})


def test_get_distances_closed_form():
    e1 = np.array([[0.0, 0.0], [1.0, 0.0]])
    e2 = np.array([[3.0, 4.0]])
    d = get_distances(e1, e2)
    assert d.shape == (2, 1)
    np.testing.assert_allclose(d[:, 0], [5.0, np.hypot(2.0, 4.0)])


def test_get_min_distances_picks_argmin():
    d = np.array([[5.0, 1.0, 9.0], [2.0, 8.0, 0.5]])
    values, indices = get_min_distances(d)
    np.testing.assert_allclose(values, [1.0, 0.5])
    assert list(indices) == [1, 2]


def test_preprocess_embeddings_shape():
    df = _emb_df(["a", "b"], [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    arr = preprocess_embeddings(df)
    assert arr.shape == (2, 3)
    np.testing.assert_allclose(arr[1], [4.0, 5.0, 6.0])


def test_min_embedding_distance_two_sets():
    generated = _emb_df(["g1", "g2"], [[0.0, 0.0], [10.0, 0.0]])
    train = _emb_df(["t_close", "t_far"], [[0.0, 1.0], [10.0, 5.0]])
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "out.csv")
        result = _min_embedding_distance(train, generated, save_path=sp)
        # g1 nearest to t_close (dist 1), g2 nearest to t_far (dist 5).
        np.testing.assert_allclose(result, [1.0, 5.0])
        out = pd.read_csv(sp).set_index("ID")
        assert out.loc["g1", "min_embedding_distance_hit"] == "t_close"
        assert out.loc["g2", "min_embedding_distance_hit"] == "t_far"
        np.testing.assert_allclose(out.loc["g1", "min_embedding_distance"], 1.0)


def test_self_mode_excludes_diagonal():
    train = _emb_df(["a", "b", "c"], [[0.0, 0.0], [1.0, 0.0], [5.0, 0.0]])
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "self.csv")
        result = _min_embedding_distance_self(train, save_path=sp)
        # a<->b dist 1, b's nearest is a (1), c's nearest is b (4). No self-hit.
        np.testing.assert_allclose(result, [1.0, 1.0, 4.0])
        out = pd.read_csv(sp).set_index("ID")
        assert out.loc["a", "min_embedding_distance_hit"] == "b"
        assert out.loc["c", "min_embedding_distance_hit"] == "b"


def test_write_topk_distances_ranking():
    distances = np.array([[9.0, 1.0, 4.0]])  # one query, three refs
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "topk.csv")
        write_topk_distances(["q"], ["r0", "r1", "r2"], distances, top_k=2, save_path=sp)
        topk = pd.read_csv(sp)
        assert list(topk.columns) == ["query_id", "rank", "neighbour_id", "score"]
        # Ascending distance: r1 (1.0) rank1, r2 (4.0) rank2.
        assert list(topk["neighbour_id"]) == ["r1", "r2"]
        assert list(topk["rank"]) == [1, 2]
        np.testing.assert_allclose(topk["score"], [1.0, 4.0])


def test_write_topk_skips_infinite_self():
    distances = np.array([[np.inf, 2.0, 3.0]])  # r0 is the excluded self
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "topk.csv")
        write_topk_distances(["q"], ["q", "r1", "r2"], distances, top_k=5, save_path=sp)
        topk = pd.read_csv(sp)
        assert "q" not in list(topk["neighbour_id"])
        assert list(topk["neighbour_id"]) == ["r1", "r2"]


def test_save_embeddings_columns():
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "s.csv")
        save_embeddings(["a", "b"], [0.1, 0.2], ["h1", "h2"], sp)
        df = pd.read_csv(sp)
        assert list(df.columns) == [
            "ID", "min_embedding_distance", "min_embedding_distance_hit"
        ]
        assert list(df["ID"]) == ["a", "b"]


def test_save_path_naming():
    assert _get_save_path("emb.csv") == "emb_min_embedding_distance.csv"
    assert _get_save_path("emb.csv", save_suffix="self") == (
        "emb_min_embedding_distance_self.csv"
    )
    assert _get_topk_save_path("emb.csv") == "emb_min_embedding_distance_topk.csv"


def test_main_train_sequences_dataframe_path():
    # main_train_sequences accepts a DataFrame directly (bypasses load_embeddings).
    train = _emb_df(["a", "b"], [[0.0, 0.0], [3.0, 0.0]])
    with tempfile.TemporaryDirectory() as d:
        sp = os.path.join(d, "self.csv")
        main_train_sequences(train, save_path=sp)
        out = pd.read_csv(sp).set_index("ID")
        np.testing.assert_allclose(out.loc["a", "min_embedding_distance"], 3.0)
        assert out.loc["a", "min_embedding_distance_hit"] == "b"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
