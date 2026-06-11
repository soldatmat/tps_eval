from __future__ import annotations

import argparse

from local_sequence_search import BACKENDS, local_sequence_search


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fast LOCAL (BLAST-style) sequence identity/similarity search of a "
        "query FASTA vs a reference, backend-pluggable over MMseqs2 (default) and "
        "DIAMOND. Writes <input>_local_sequence_search.csv (keyed by ID) with best-hit "
        "local identity/similarity/coverage. With --top_k N, also writes "
        "<input>_local_sequence_search_topk.csv (columns query_id,rank,neighbour_id,score; "
        "score = identity PERCENT in [0,100], LARGER = closer). Complements the GLOBAL "
        "max_sequence_identity tool; does not replace it."
    )
    parser.add_argument("fasta_path", help="Query FASTA file.")
    parser.add_argument(
        "train_path",
        nargs="?",
        default=None,
        help="Optional reference FASTA. If omitted, self mode "
        "(each query's best hit / neighbours exclude itself).",
    )
    parser.add_argument(
        "--backend",
        choices=BACKENDS,
        default="mmseqs2",
        help="Search backend (default: mmseqs2).",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=None,
        help="If >= 1, also emit the top-k nearest reference neighbours per query.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Backend threads (default 4).")
    parser.add_argument(
        "--sensitivity",
        default=None,
        help="Backend sensitivity knob (mmseqs2 -s value e.g. 7.5; diamond flag name "
        "e.g. very-sensitive). Default: backend default.",
    )
    parser.add_argument("--save_path", default=None, help="Metric CSV path.")
    parser.add_argument("--topk_save_path", default=None, help="Top-k CSV path.")
    args = parser.parse_args()

    local_sequence_search(
        args.fasta_path,
        train_path=args.train_path,
        backend=args.backend,
        top_k=args.top_k,
        threads=args.threads,
        sensitivity=args.sensitivity,
        save_path=args.save_path,
        topk_save_path=args.topk_save_path,
    )


if __name__ == "__main__":
    main()
