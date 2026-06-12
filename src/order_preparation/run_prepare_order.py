"""argv entry point for order preparation. See ``prepare_order.py`` for the logic.

Examples
--------
    # Batch from a FASTA (default: yeast, Type 3 overhangs):
    python run_prepare_order.py designs.fasta

    # From a CSV, custom output prefix and organism:
    python run_prepare_order.py designs.csv -o my_order --organism yeast

    # A single inline sequence:
    python run_prepare_order.py --sequence MGRSY...PIPL --id design1
"""
from __future__ import annotations

import argparse
import sys

from codon_optimization import (
    DEFAULT_GC_MAX,
    DEFAULT_GC_MIN,
    DEFAULT_GC_WINDOW,
    DEFAULT_MAX_HOMOPOLYMER,
    DEFAULT_METHOD,
    DEFAULT_ORGANISM,
    DEFAULT_SEED,
)
from overhangs import DEFAULT_OVERHANG, OVERHANGS
from prepare_order import prepare_one, prepare_order


def main() -> int:
    p = argparse.ArgumentParser(
        description="Codon-optimize protein designs and add Golden Gate overhangs for ordering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_path", nargs="?", help="FASTA or CSV/TSV of protein designs.")
    p.add_argument("--sequence", help="A single amino-acid sequence (instead of input_path).")
    p.add_argument("--id", default="design", help="ID for --sequence mode.")
    p.add_argument("-o", "--output_prefix", help="Output path prefix (default: input without suffix).")
    p.add_argument("--organism", default=DEFAULT_ORGANISM, help="Target organism for codon usage.")
    p.add_argument(
        "--overhang_type", default=DEFAULT_OVERHANG, choices=list(OVERHANGS),
        metavar="TYPE", help="Golden Gate overhang type (see Plasmid_Generator.xlsx).",
    )
    p.add_argument("--id-column", dest="id_column", help="ID column name (CSV input).")
    p.add_argument("--seq-column", dest="seq_column", help="Amino-acid column name (CSV input).")
    p.add_argument(
        "--method", default=DEFAULT_METHOD,
        help="DNAChisel CodonOptimize method (match_codon_usage | use_best_codon | harmonize_rca).",
    )
    p.add_argument(
        "--max_homopolymer", type=int, default=DEFAULT_MAX_HOMOPOLYMER,
        help="Longest allowed single-nucleotide run (0 disables the cap).",
    )
    p.add_argument("--gc_min", type=float, default=DEFAULT_GC_MIN, help="GC-window lower bound (fraction).")
    p.add_argument("--gc_max", type=float, default=DEFAULT_GC_MAX, help="GC-window upper bound (fraction).")
    p.add_argument(
        "--gc_window", type=int, default=DEFAULT_GC_WINDOW,
        help="GC sliding-window size in bp (0 disables the GC window).",
    )
    p.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help="RNG seed for reproducible codon sampling (use -1 for nondeterministic).",
    )
    args = p.parse_args()
    seed = None if args.seed is not None and args.seed < 0 else args.seed

    quality = dict(
        method=args.method, max_homopolymer=args.max_homopolymer,
        gc_min=args.gc_min, gc_max=args.gc_max, gc_window=args.gc_window, seed=seed,
    )

    if args.sequence:
        row = prepare_one(args.sequence, args.organism, args.overhang_type, **quality)
        if row["warnings"]:
            print(f"[WARN] {args.id}: {row['warnings']}", file=sys.stderr)
        print(f"{args.id},{row['ordered_sequence']}")
        return 0

    if not args.input_path:
        p.error("provide an input_path or --sequence")

    prepare_order(
        args.input_path,
        output_prefix=args.output_prefix,
        organism=args.organism,
        overhang_type=args.overhang_type,
        id_column=args.id_column,
        seq_column=args.seq_column,
        **quality,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
