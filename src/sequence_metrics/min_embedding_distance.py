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

    return min_dist.tolist()


def _min_embedding_distance_self(
    train_df: pd.DataFrame,
    *,
    save_path: Optional[str] = None,
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

    return min_dist.tolist()


def _get_save_path(embeddings_path: str, *, save_suffix: Optional[str] = None) -> str:
    extension = embeddings_path.split(".")[-1]
    base_path = embeddings_path[: -len(extension) - 1]
    suffix = "" if save_suffix is None else f"_{save_suffix}"
    return f"{base_path}_min_embedding_distance{suffix}.csv"


def min_embedding_distance(
    embeddings_path: str,
    *,
    train_embeddings_path: Optional[str] = None,
    save: bool = True,
):
    """Compute per-sequence minimum embedding distance.

    Args:
        embeddings_path: CSV with evaluated embeddings.
        train_embeddings_path: Optional CSV with reference embeddings.
        save: Save output CSV when True.
    """
    if train_embeddings_path is None:
        main_train_sequences(embeddings_path, save=save)
    else:
        main_generated_sequences(embeddings_path, train_embeddings_path, save=save)


def main_generated_sequences(
    generated_embeddings_path: str,
    train_embeddings_path: str,
    *,
    save: bool = True,
):
    train_df = load_embeddings(train_embeddings_path)
    generated_df = load_embeddings(generated_embeddings_path)

    save_path = _get_save_path(generated_embeddings_path)
    _min_embedding_distance(
        train_df,
        generated_df,
        save_path=save_path if save else None,
    )


def main_train_sequences(
    embeddings_path,
    *,
    save: bool = True,
    save_path: Optional[str] = None,
):
    if isinstance(embeddings_path, pd.DataFrame):
        train_df = embeddings_path
        _min_embedding_distance_self(train_df, save_path=save_path)
        return

    train_df = load_embeddings(embeddings_path)
    resolved_save_path = _get_save_path(embeddings_path, save_suffix="self")
    _min_embedding_distance_self(
        train_df,
        save_path=resolved_save_path if save else None,
    )
