from __future__ import annotations

import sys

from min_embedding_distance import min_embedding_distance


def main() -> None:
    num_args = len(sys.argv) - 1
    if num_args == 1:
        embeddings_path = sys.argv[1]
        min_embedding_distance(embeddings_path, save=True)
    elif num_args == 2:
        embeddings_path = sys.argv[1]
        train_embeddings_path = sys.argv[2]
        min_embedding_distance(
            embeddings_path,
            train_embeddings_path=train_embeddings_path,
            save=True,
        )
    else:
        raise ValueError(f"Invalid number of arguments. Expected 1 or 2, got {num_args}.")


if __name__ == "__main__":
    main()
