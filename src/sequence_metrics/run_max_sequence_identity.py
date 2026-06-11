from __future__ import annotations

import argparse

from max_sequence_identity import max_sequence_identity


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-query maximum sequence identity to a reference set. "
        "With --top_k N, additionally writes <input>_max_sequence_identity_topk.csv "
        "(columns query_id,rank,neighbour_id,score) where score is identity PERCENT "
        "in [0, 100] (LARGER = closer). Default single-best output is unchanged."
    )
    parser.add_argument("fasta_path", help="Path to the FASTA file to evaluate.")
    parser.add_argument(
        "train_path",
        nargs="?",
        default=None,
        help="Optional reference FASTA. If omitted, self mode "
        "(each query's neighbours exclude itself).",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=None,
        help="If >= 1, also emit the top-k nearest reference hits per query.",
    )
    args = parser.parse_args()

    max_sequence_identity(
        args.fasta_path,
        train_path=args.train_path,
        top_k=args.top_k,
    )


if __name__ == "__main__":
    main()
