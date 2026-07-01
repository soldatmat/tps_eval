#!/usr/bin/env python3
"""argv entry for the composite `select` selection driver.

Usage:
    # from an already-merged wide table:
    python run_select.py --merged merged.csv --spec select_phase3.json \
        --output_prefix phase3 [--fasta seed.fasta]

    # or merge the per-tool CSVs inline, then select:
    python run_select.py --entries <csv|dir|glob> [...] --spec select_phase1.json \
        --output_prefix phase1 --fasta seed.fasta

Writes <prefix>_survivors.csv, <prefix>_survivors.fasta, and <prefix>_manifest.md.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from io_fasta import read_fasta_map
from merge import merge_metrics
from select_designs import select_and_write


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--merged", help="Pre-merged wide metric CSV (keyed by ID).")
    src.add_argument("--entries", nargs="+",
                     help="Per-tool CSVs / dirs / globs to merge inline before selecting.")
    p.add_argument("--spec", required=True, help="Selection spec JSON.")
    p.add_argument("--output_prefix", required=True, help="Output path prefix.")
    p.add_argument("--fasta", default=None,
                   help="Seed FASTA to source sequences from when the merged table lacks a "
                        "'sequence' column (needed to emit the survivor FASTA).")
    p.add_argument("--title", default="Selection", help="Manifest title.")
    args = p.parse_args()

    with open(args.spec) as fh:
        spec = json.load(fh)
    df = pd.read_csv(args.merged) if args.merged else merge_metrics(args.entries)
    fasta_map = read_fasta_map(args.fasta) if args.fasta else None
    select_and_write(df, spec, args.output_prefix, fasta_map=fasta_map, title=args.title)
    return 0


if __name__ == "__main__":
    sys.exit(main())
