from __future__ import annotations

import argparse

from foldseek_swissprot_search import foldseek_swissprot_search


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Broad structure homology search of designs vs AlphaFold-Swiss-Prot "
        "(foldseek easy-search), reporting the top hit across all proteins and whether "
        "the top-N hits are terpene synthases. Writes a CSV keyed by ID."
    )
    parser.add_argument("structs_dir", help="Directory of generated structures "
                        "(af_output dir or flat dir of .pdb/.cif).")
    parser.add_argument("afdb_db", help="foldseek AlphaFold/Swiss-Prot DB path.")
    parser.add_argument("tps_accessions_path", help="Committed TPS accession list "
                        "(data/reference/tps_uniprot_accessions.txt).")
    parser.add_argument("--save_path", default=None,
                        help="Output CSV path (default: <structs_dir>_foldseek_swissprot_search.csv).")
    parser.add_argument("--top_n", type=int, default=25, help="Top-N hits per query (default 25).")
    parser.add_argument("--max_seqs", type=int, default=300,
                        help="foldseek --max-seqs prefilter (default 300).")
    args = parser.parse_args()

    foldseek_swissprot_search(
        args.structs_dir,
        args.afdb_db,
        args.tps_accessions_path,
        save_path=args.save_path,
        top_n=args.top_n,
        max_seqs=args.max_seqs,
    )


if __name__ == "__main__":
    main()
