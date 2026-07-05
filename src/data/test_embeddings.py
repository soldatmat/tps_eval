from __future__ import annotations

"""Self-contained tests for data/embeddings.py (embedding CSV loading).

Run from this directory:
    cd src/data && python test_embeddings.py
or under pytest:
    cd src/data && python -m pytest test_embeddings.py -q

Synthetic in-memory CSVs only. No torch is required (this loader is pure pandas +
ast). Covers the cell parser, the id-column auto-detection, both the ``embedding``
list-column layout and the wide feature-column layout, and a CSV round-trip.
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.embeddings import _parse_embedding_cell, load_embeddings  # noqa: E402


def _tmp_csv(df: pd.DataFrame) -> str:
    tmp = tempfile.mkdtemp(prefix="embeddings_")
    path = os.path.join(tmp, "emb.csv")
    df.to_csv(path, index=False)
    return path


def test_parse_embedding_cell():
    assert _parse_embedding_cell([1, 2, 3]) == [1.0, 2.0, 3.0]
    assert _parse_embedding_cell("[1.0, 2.5]") == [1.0, 2.5]
    try:
        _parse_embedding_cell(3.14)   # unsupported scalar type
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for unsupported cell type")
    print("ok parse_embedding_cell")


def test_load_embeddings_list_column_roundtrip():
    df = pd.DataFrame({"ID": ["a", "b"], "embedding": [[1.0, 2.0], [3.0, 4.0]]})
    path = _tmp_csv(df)   # embedding serialized as "[1.0, 2.0]" strings
    out = load_embeddings(path)
    assert list(out.columns) == ["ID", "embedding"]
    assert list(out["ID"]) == ["a", "b"]
    assert out["embedding"].tolist() == [[1.0, 2.0], [3.0, 4.0]]
    print("ok load_embeddings_list_column_roundtrip")


def test_load_embeddings_wide_feature_columns():
    df = pd.DataFrame({"id": ["x", "y"], "f0": [1.0, 3.0], "f1": [2.0, 4.0]})
    out = load_embeddings(_tmp_csv(df))
    assert list(out["ID"]) == ["x", "y"]
    assert out["embedding"].tolist() == [[1.0, 2.0], [3.0, 4.0]]
    print("ok load_embeddings_wide_feature_columns")


def test_id_column_is_stringified():
    df = pd.DataFrame({"id": [1, 2], "embedding": [[0.1], [0.2]]})
    out = load_embeddings(_tmp_csv(df))
    assert list(out["ID"]) == ["1", "2"]   # ids coerced to str
    print("ok id_column_is_stringified")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
