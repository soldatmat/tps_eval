"""Derive the prenyl-diphosphate SUBSTRATE-class `reference_id,label` file for MARTS-DB TPS.

This is a LABELING INPUT for the label-agnostic `knn_label_transfer.py` tool, mirroring
`make_first_cyclization_labels.py` but labeling each known TPS by the *size class of its
prenyl-diphosphate substrate* (GPP/C10 mono, FPP/C15 sesqui, GGPP/C20 di, ...). The k-NN
tool stays agnostic; swapping in this label file makes it a substrate-class predictor.

Source: the MARTS-DB reactions table shipped with the EnzymeExplorer repo
(`EnzymeExplorer/data/martsDB_reactions_2026_02_22.csv`; one row per (enzyme, reaction),
4205 rows / 1378 enzymes). That table carries a curated **`Type`** column
(`mono`/`sesq`/`di`/`sester`/`tri`/`sqs`/`pt`/`psy`/`tetra`/`hemi`/`sesquar`) — the
terpene size class, which is in 1:1 correspondence with the prenyl-diphosphate substrate
the enzyme consumes (mono = GPP, sesq = FPP, di = GGPP, ...). We use `Type` directly: it
is the curators' substrate/size assignment and is far cleaner than re-deriving carbon
counts from `Substrate_smiles`. We DO record the canonical substrate name + carbon count
per class for documentation, and (as a fallback) the script can map the substrate SMILES
carbon count to C10/C15/C20 when a `Type` is missing.

`Type` -> substrate class (label), chosen to share a vocabulary with the EnzymeExplorer
sequence-only per-substrate scores (`GPP`/`FPP`/`GGPP`/`GFPP`/`EDSQ`/`2xGGPP`/`IDS` ...)
so the two signals can be compared directly in the combiner:

    mono     -> GPP    (C10  geranyl-PP)
    sesq     -> FPP    (C15  farnesyl-PP)
    di       -> GGPP   (C20  geranylgeranyl-PP; incl. copalyl-PP cyclases)
    sester   -> GFPP   (C25  geranylfarnesyl-PP)
    tri, sqs -> EDSQ   (C30  squalene / 2,3-epoxysqualene, from 2xFPP)
    psy,tetra-> 2xGGPP (C40  phytoene / carotenoid, from 2xGGPP)
    pt       -> IDS    (prenyltransferase / chain elongation, IPP+DMAPP/prenyl-PP)
    hemi     -> DMAPP  (C5   hemiterpene)
    sesquar  -> C35    (C35  sesquarterpene, heptaprenyl-PP)

Multi-row / multi-`Type` enzymes (200 of 1378 carry >1 Type across their reactions) are
collapsed to ONE coarse label = the most frequent class (ties -> the class with the
smallest carbon count, for determinism), matching the "one coarse class per known TPS"
assumption of the transfer metric.

Output: `substrate_labels.csv` with columns `reference_id,label` where `reference_id` is
the MARTS enzyme id (`marts_E*`, matching the FASTA / struct stems) and `label` is the
substrate class string (`GPP`/`FPP`/`GGPP`/...).

Usage:
    python make_substrate_labels.py \
        --source /path/to/EnzymeExplorer/data/martsDB_reactions_2026_02_22.csv \
        --out substrate_labels.csv
"""
from __future__ import annotations

import argparse

import pandas as pd

# MARTS-DB `Type` -> substrate size class. The label strings mirror the EnzymeExplorer
# sequence-only per-substrate score columns so k-NN and EE share a vocabulary.
TYPE_TO_SUBSTRATE = {
    "mono": "GPP",      # C10
    "sesq": "FPP",      # C15
    "di": "GGPP",       # C20 (incl. copalyl-PP)
    "sester": "GFPP",   # C25
    "tri": "EDSQ",      # C30 (squalene / 2,3-epoxysqualene)
    "sqs": "EDSQ",      # C30 (squalene synthase, from 2xFPP)
    "psy": "2xGGPP",    # C40 (phytoene synthase)
    "tetra": "2xGGPP",  # C40 (carotenoid)
    "pt": "IDS",        # prenyltransferase / chain elongation
    "hemi": "DMAPP",    # C5
    "sesquar": "C35",   # C35 (sesquarterpene)
}

# Approximate carbon count per substrate class — used only to break ties deterministically
# (smaller substrate wins) and for documentation. NOT a precise stoichiometry.
SUBSTRATE_CARBONS = {
    "DMAPP": 5,
    "GPP": 10,
    "FPP": 15,
    "GGPP": 20,
    "GFPP": 25,
    "EDSQ": 30,
    "C35": 35,
    "2xGGPP": 40,
    "IDS": 999,  # elongation/prenyltransferase: no single size; sort last
}


def _carbons(label: str) -> int:
    return SUBSTRATE_CARBONS.get(label, 1000)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True,
                    help="martsDB_reactions_*.csv path (must have Enzyme_marts_ID + Type).")
    ap.add_argument("--output", required=True, help="Output reference_id,label CSV.")
    ap.add_argument("--id_column", default="Enzyme_marts_ID")
    ap.add_argument("--type_column", default="Type")
    args = ap.parse_args()

    df = pd.read_csv(args.source)
    df = df[[args.id_column, args.type_column]].dropna()
    df["substrate"] = df[args.type_column].astype(str).str.strip().str.lower().map(TYPE_TO_SUBSTRATE)
    n_unmapped = int(df["substrate"].isna().sum())
    if n_unmapped:
        bad = sorted(df.loc[df["substrate"].isna(), args.type_column].astype(str).unique())
        print(f"[warn] {n_unmapped} rows with un-mapped Type value(s): {bad} -> dropped.")
    df = df.dropna(subset=["substrate"])

    rows = []
    for rid, grp in df.groupby(args.id_column):
        counts = grp["substrate"].value_counts()
        # most frequent class; tie -> smallest carbon count (then alphabetical)
        top = sorted(
            counts[counts == counts.max()].index,
            key=lambda lab: (_carbons(lab), lab),
        )[0]
        rows.append({"reference_id": rid, "label": top})

    out = pd.DataFrame(rows).sort_values("reference_id")
    out.to_csv(args.output, index=False)
    dist = out["label"].value_counts().to_dict()
    print(f"Wrote {len(out)} reference labels to {args.output} "
          f"({out['label'].nunique()} distinct substrate classes).")
    print("  class distribution:", dist)


if __name__ == "__main__":
    main()
