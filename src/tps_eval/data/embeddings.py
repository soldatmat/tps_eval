from __future__ import annotations

import ast
from typing import List

import pandas as pd


def _parse_embedding_cell(value) -> List[float]:
    if isinstance(value, list):
        return [float(v) for v in value]
    if isinstance(value, str):
        parsed = ast.literal_eval(value)
        return [float(v) for v in parsed]
    raise ValueError(f"Unsupported embedding value type: {type(value)!r}")


def load_embeddings(file_path: str) -> pd.DataFrame:
    """Load embedding vectors into a DataFrame with columns ['ID', 'embedding']."""
    original_df = pd.read_csv(file_path)

    if "id" in original_df.columns:
        id_col = "id"
    elif "ID" in original_df.columns:
        id_col = "ID"
    else:
        id_col = original_df.columns[0]

    if "embedding" in original_df.columns:
        sequence_embeddings = [_parse_embedding_cell(v) for v in original_df["embedding"].tolist()]
    else:
        feature_df = original_df.drop(columns=[id_col])
        sequence_embeddings = []
        for _, row in feature_df.iterrows():
            sequence_embeddings.append([float(v) for v in row.tolist()])

    embedding_df = pd.DataFrame(
        {
            "ID": original_df[id_col].astype(str),
            "embedding": sequence_embeddings,
        }
    )
    return embedding_df
