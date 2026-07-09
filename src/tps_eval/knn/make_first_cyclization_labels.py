"""Derive the first-cyclization `reference_id,label` file for MARTS-DB TPS.

This is the LABELING INPUT for `knn_label_transfer.py` under today's labeling
(first-cyclization class, which aligns with the generation-side conditioning). The
tool itself is label-agnostic; this script just produces one valid label file.

Source: the companion `tps-first-cyclization-knn` project's curated table
`data/TPS_first_cyclization.csv` (one row per (enzyme, product); 1349 rows / 991
enzymes; `First_cyclization_product_id` in 0..21). That table is the same labeled
MARTS-DB set the kNN-conditioning-fidelity work used. A few enzymes carry more than
one first-cyclization class across their products (multi-product enzymes); we collapse
to ONE coarse label per enzyme = the most frequent class (ties -> smallest class id,
for determinism), matching the "one coarse class per known TPS" assumption of the
transfer metric.

Output: `first_cyclization_labels.csv` with columns `reference_id,label` where
`reference_id` is the MARTS enzyme id (`marts_E*`, matching the FASTA / struct stems)
and `label` is the integer first-cyclization class id.

Usage:
    python make_first_cyclization_labels.py \
        --source /path/to/tps-first-cyclization-knn/data/TPS_first_cyclization.csv \
        --out first_cyclization_labels.csv
"""
from __future__ import annotations

import argparse

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, help="TPS_first_cyclization.csv path.")
    ap.add_argument("--output", required=True, help="Output reference_id,label CSV.")
    ap.add_argument("--id_column", default="Enzyme_marts_ID")
    ap.add_argument("--label_column", default="First_cyclization_product_id")
    args = ap.parse_args()

    df = pd.read_csv(args.source)
    df = df[[args.id_column, args.label_column]].dropna()
    df[args.label_column] = df[args.label_column].astype(int)

    rows = []
    for rid, grp in df.groupby(args.id_column):
        # most frequent class; tie -> smallest class id
        counts = grp[args.label_column].value_counts()
        top = counts[counts == counts.max()].index.min()
        rows.append({"reference_id": rid, "label": int(top)})

    out = pd.DataFrame(rows).sort_values("reference_id")
    out.to_csv(args.output, index=False)
    print(f"Wrote {len(out)} reference labels to {args.output} "
          f"({out['label'].nunique()} distinct classes).")


if __name__ == "__main__":
    main()
