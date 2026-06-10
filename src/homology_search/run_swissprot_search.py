from __future__ import annotations

import argparse

from swissprot_search import swissprot_search


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Broad sequence homology search of designs vs Swiss-Prot (DIAMOND "
        "blastp), reporting the top hit across all proteins and whether the top-N hits "
        "are terpene synthases. Writes a CSV keyed by ID."
    )
    parser.add_argument("fasta_path", help="Design FASTA file.")
    parser.add_argument("diamond_db", help="DIAMOND DB built from uniprot_sprot.fasta.")
    parser.add_argument("tps_accessions_path", help="Committed TPS accession list "
                        "(src/homology_search/tps_uniprot_accessions.txt).")
    parser.add_argument("--save_path", default=None,
                        help="Output CSV path (default: <fasta>_swissprot_search.csv).")
    parser.add_argument("--top_n", type=int, default=25, help="Top-N hits per query (default 25).")
    parser.add_argument("--threads", type=int, default=4, help="DIAMOND threads (default 4).")
    parser.add_argument("--sensitivity", default="very-sensitive",
                        help="DIAMOND sensitivity flag, e.g. sensitive, very-sensitive, "
                        "ultra-sensitive (default very-sensitive).")
    args = parser.parse_args()

    swissprot_search(
        args.fasta_path,
        args.diamond_db,
        args.tps_accessions_path,
        save_path=args.save_path,
        top_n=args.top_n,
        threads=args.threads,
        sensitivity=args.sensitivity,
    )


if __name__ == "__main__":
    main()
