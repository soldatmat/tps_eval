#!/usr/bin/env python3
"""argv entry for the merge_metrics selection primitive.

Usage:
    python run_merge.py --entries <csv|dir|glob> [<csv|dir|glob> ...] --output merged.csv

Each entry is a per-tool CSV, a directory of them, or a glob. They are merged into one
wide table keyed by ID (see merge.py for the conventions).
"""
from __future__ import annotations

import argparse
import sys

from merge import merge_metrics, write_merged


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--entries", nargs="+", required=True,
                   help="Per-tool CSVs / directories / globs to merge (keyed by ID).")
    p.add_argument("--output", required=True, help="Output merged CSV path.")
    args = p.parse_args()
    df = merge_metrics(args.entries)
    write_merged(df, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
