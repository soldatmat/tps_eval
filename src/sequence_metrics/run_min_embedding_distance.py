from __future__ import annotations

import argparse

from min_embedding_distance import min_embedding_distance


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-query minimum ESM-embedding distance to a reference set. "
        "With --top_k N, additionally writes <input>_min_embedding_distance_topk.csv "
        "(columns query_id,rank,neighbour_id,score) where score is the ESM-embedding "
        "L2 distance (SMALLER = closer). Default single-best output is unchanged."
    )
    parser.add_argument("embeddings_path", help="CSV with the embeddings to evaluate.")
    parser.add_argument(
        "train_embeddings_path",
        nargs="?",
        default=None,
        help="Optional reference embeddings CSV. If omitted, self mode "
        "(each query's neighbours exclude itself).",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=None,
        help="If >= 1, also emit the top-k nearest reference neighbours per query.",
    )
    args = parser.parse_args()

    min_embedding_distance(
        args.embeddings_path,
        train_embeddings_path=args.train_embeddings_path,
        save=True,
        top_k=args.top_k,
    )


if __name__ == "__main__":
    main()
