from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.embeddings import load_embeddings


def preprocess_embeddings(embeddings_df: pd.DataFrame) -> np.ndarray:
    return np.array(embeddings_df["embedding"].tolist(), dtype=float)


def get_distances(embeddings1: np.ndarray, embeddings2: np.ndarray) -> np.ndarray:
    diffs = embeddings1[:, None, :] - embeddings2[None, :, :]
    return np.linalg.norm(diffs, axis=2)


def get_min_distances(distances: np.ndarray):
    indices = np.argmin(distances, axis=1)
    values = distances[np.arange(distances.shape[0]), indices]
    return values, indices


def write_topk_distances(
    query_ids,
    train_ids,
    distances: np.ndarray,
    top_k: int,
    save_path: str,
) -> None:
    """Write the top-k nearest reference neighbours (SMALLEST distance) per query.

    Tidy CSV with columns query_id,rank,neighbour_id,score. ``score`` is the
    ESM-embedding L2 distance (SMALLER = closer). Self-exclusion is expected to
    be already applied to ``distances`` (e.g. diagonal set to inf in self mode).
    """
    k = min(top_k, distances.shape[1])
    train_ids = list(train_ids)
    rows = []
    for q in range(distances.shape[0]):
        row = distances[q]
        # Ascending by distance; np.argsort is stable so ties break on train index.
        order = np.argsort(row, kind="stable")
        rank = 0
        for j in order:
            score = row[j]
            if not np.isfinite(score):
                continue  # excluded (e.g. self in self mode)
            rank += 1
            rows.append(
                {
                    "query_id": str(query_ids[q]).split()[0] if str(query_ids[q]) else query_ids[q],
                    "rank": rank,
                    "neighbour_id": str(train_ids[j]).split()[0],
                    "score": float(score),
                }
            )
            if rank >= k:
                break
    pd.DataFrame(rows, columns=["query_id", "rank", "neighbour_id", "score"]).to_csv(
        save_path, index=False
    )


def save_embeddings(
    ids,
    min_dist,
    min_dist_hits,
    save_path: str,
) -> None:
    df = pd.DataFrame(
        {
            "ID": ids,
            "min_embedding_distance": min_dist,
            "min_embedding_distance_hit": min_dist_hits,
        }
    )
    df.to_csv(save_path, index=False)


def _min_embedding_distance(
    train_df: pd.DataFrame,
    generated_df: pd.DataFrame,
    *,
    save_path: Optional[str] = None,
    top_k: Optional[int] = None,
    topk_save_path: Optional[str] = None,
):
    generated_embeddings = preprocess_embeddings(generated_df)
    train_embeddings = preprocess_embeddings(train_df)

    distances = get_distances(generated_embeddings, train_embeddings)
    min_dist, min_dist_index = get_min_distances(distances)
    min_dist_hits = train_df.iloc[min_dist_index]["ID"].tolist()

    if save_path is not None:
        save_embeddings(
            generated_df["ID"].tolist(),
            min_dist.tolist(),
            min_dist_hits,
            save_path,
        )

    if top_k is not None and top_k >= 1 and topk_save_path is not None:
        write_topk_distances(
            generated_df["ID"].tolist(),
            train_df["ID"].tolist(),
            distances,
            top_k,
            topk_save_path,
        )

    return min_dist.tolist()


def _min_embedding_distance_self(
    train_df: pd.DataFrame,
    *,
    save_path: Optional[str] = None,
    top_k: Optional[int] = None,
    topk_save_path: Optional[str] = None,
):
    train_embeddings = preprocess_embeddings(train_df)
    distances = get_distances(train_embeddings, train_embeddings)
    np.fill_diagonal(distances, np.inf)

    min_dist, min_dist_index = get_min_distances(distances)
    min_dist_hits = train_df.iloc[min_dist_index]["ID"].tolist()

    if save_path is not None:
        save_embeddings(
            train_df["ID"].tolist(),
            min_dist.tolist(),
            min_dist_hits,
            save_path,
        )

    if top_k is not None and top_k >= 1 and topk_save_path is not None:
        # Diagonal is inf -> self excluded from each query's neighbour list.
        write_topk_distances(
            train_df["ID"].tolist(),
            train_df["ID"].tolist(),
            distances,
            top_k,
            topk_save_path,
        )

    return min_dist.tolist()


def _get_save_path(embeddings_path: str, *, save_suffix: Optional[str] = None) -> str:
    extension = embeddings_path.split(".")[-1]
    base_path = embeddings_path[: -len(extension) - 1]
    suffix = "" if save_suffix is None else f"_{save_suffix}"
    return f"{base_path}_min_embedding_distance{suffix}.csv"


def _get_topk_save_path(embeddings_path: str) -> str:
    extension = embeddings_path.split(".")[-1]
    base_path = embeddings_path[: -len(extension) - 1]
    return f"{base_path}_min_embedding_distance_topk.csv"


def min_embedding_distance(
    embeddings_path: str,
    *,
    train_embeddings_path: Optional[str] = None,
    save: bool = True,
    top_k: Optional[int] = None,
):
    """Compute per-sequence minimum embedding distance.

    Args:
        embeddings_path: CSV with evaluated embeddings.
        train_embeddings_path: Optional CSV with reference embeddings.
        save: Save output CSV when True.
        top_k: When >= 1, also write <input>_min_embedding_distance_topk.csv
            (columns query_id,rank,neighbour_id,score). ``score`` is the ESM-
            embedding L2 distance (SMALLER = closer). In self mode each query
            excludes itself.
    """
    if train_embeddings_path is None:
        main_train_sequences(embeddings_path, save=save, top_k=top_k)
    else:
        main_generated_sequences(
            embeddings_path, train_embeddings_path, save=save, top_k=top_k
        )


def main_generated_sequences(
    generated_embeddings_path: str,
    train_embeddings_path: str,
    *,
    save: bool = True,
    top_k: Optional[int] = None,
):
    train_df = load_embeddings(train_embeddings_path)
    generated_df = load_embeddings(generated_embeddings_path)

    save_path = _get_save_path(generated_embeddings_path)
    _min_embedding_distance(
        train_df,
        generated_df,
        save_path=save_path if save else None,
        top_k=top_k,
        topk_save_path=_get_topk_save_path(generated_embeddings_path),
    )


def main_train_sequences(
    embeddings_path,
    *,
    save: bool = True,
    save_path: Optional[str] = None,
    top_k: Optional[int] = None,
    topk_save_path: Optional[str] = None,
):
    if isinstance(embeddings_path, pd.DataFrame):
        train_df = embeddings_path
        _min_embedding_distance_self(
            train_df,
            save_path=save_path,
            top_k=top_k,
            topk_save_path=topk_save_path,
        )
        return

    train_df = load_embeddings(embeddings_path)
    resolved_save_path = _get_save_path(embeddings_path, save_suffix="self")
    _min_embedding_distance_self(
        train_df,
        save_path=resolved_save_path if save else None,
        top_k=top_k,
        topk_save_path=_get_topk_save_path(embeddings_path),
    )
