#!/usr/bin/env python3
"""Render a first-cyclization-class-coloured 2D landscape map of a protein set.

Standalone visualization runner (not part of the per-design eval orchestrator).
Consumes a representation produced by one of the per-domain packages and lays it
out with one or more DR backends, one panel per method.

Two input modes:
  --features CSV   : feature matrix, first col `id` (Enzyme_marts_ID) then dims.
                     Methods: pca, tsne, umap, pacmap.
  --pairs TSV      : all-vs-all similarity table (+ --pair-cols, --sim-col);
                     distance = 1 - sim, missing pairs -> 1. Methods: umap, tsne, pcoa
                     (all precomputed-distance).

Class labels (for colouring) via one of:
  --label-col COL              : label already in the features CSV.
  --labels-parallel CSV --class-col COL : labels by ROW POSITION (CSV row-aligned
                                          with the features CSV; for the 1349-row
                                          ESM matrix).
  --labels-join CSV --id-col COL --class-col COL : join by id (first class per id;
                                          for unique-id sets / --pairs).

Example:
  python scripts/run_visualization.py --features data/visualization/active_site_features.csv \
    --label-col product_class_id --methods umap --title "Active-site features" \
    --out data/visualization/active_site_umap.png
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.visualization import dimensionality_reduction as dr  # noqa: E402
from src.visualization import landscape_map  # noqa: E402


def build_distance(df, ids, qcol, tcol, scol):
    pos = {i: k for k, i in enumerate(ids)}
    n = len(ids)
    S = np.zeros((n, n), float)
    for qi, ti, si in zip(df[qcol].astype(str), df[tcol].astype(str),
                          df[scol].astype(float)):
        i, j = pos.get(qi), pos.get(ti)
        if i is None or j is None:
            continue
        if si > S[i, j]:
            S[i, j] = si
            S[j, i] = si
    S = np.maximum(S, S.T)
    np.fill_diagonal(S, 1.0)
    D = np.clip(1.0 - S, 0.0, 1.0)
    np.fill_diagonal(D, 0.0)
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features")
    ap.add_argument("--pairs")
    ap.add_argument("--pair-cols", help="comma list of column names for --pairs (no header)")
    ap.add_argument("--sim-col")
    ap.add_argument("--label-col")
    ap.add_argument("--labels-parallel")
    ap.add_argument("--labels-join")
    ap.add_argument("--id-col", default="Enzyme_marts_ID")
    ap.add_argument("--class-col", default="First_cyclization_product_id")
    ap.add_argument("--exclude-cols", default="", help="comma list of extra (non-feature) numeric columns to drop")
    ap.add_argument("--methods", default="umap")
    ap.add_argument("--title", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--footnote", default="")
    args = ap.parse_args()
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    panels = []
    if args.features:
        df = pd.read_csv(args.features)
        ids = df["id"].astype(str).to_numpy() if "id" in df.columns else df.iloc[:, 0].astype(str).to_numpy()
        # numeric feature columns only
        meta = {"id", args.id_col, args.class_col, args.label_col}
        meta |= {c.strip() for c in args.exclude_cols.split(",") if c.strip()}
        feat_cols = [c for c in df.columns if c not in meta and pd.api.types.is_numeric_dtype(df[c])]
        X = df[feat_cols].to_numpy(dtype=float)
        keep = ~np.isnan(X).any(axis=1)
        X, ids = X[keep], ids[keep]
        y = _labels_for_features(df[keep], ids, args)
        print(f"features: {X.shape} | {len(feat_cols)} dims | n={len(ids)}")
        for m in methods:
            if m == "pca":
                c, (p1, p2) = dr.pca_2d(X); panels.append((c, f"PCA  ({p1:.1f}%, {p2:.1f}%)"))
            elif m == "tsne":
                panels.append((dr.tsne_2d(X), "t-SNE  (perplexity 30)"))
            elif m == "umap":
                panels.append((dr.umap_2d(X), "UMAP"))
            elif m == "pacmap":
                panels.append((dr.pacmap_2d(X), "PaCMAP"))
            else:
                raise SystemExit(f"method {m} not valid for --features")
    elif args.pairs:
        cols = [c.strip() for c in args.pair_cols.split(",")]
        df = pd.read_csv(args.pairs, sep="\t", header=None, names=cols)
        ids = sorted(set(df[cols[0]].astype(str)) | set(df[cols[1]].astype(str)))
        D = build_distance(df, ids, cols[0], cols[1], args.sim_col)
        print(f"pairs: {len(df)} | ids={len(ids)} | mean off-diag dist "
              f"{D[~np.eye(len(ids),dtype=bool)].mean():.3f}")
        y = _labels_join(np.array(ids), args.labels_join or args.labels_parallel, args)
        for m in methods:
            if m == "umap":
                panels.append((dr.umap_2d(D, precomputed=True), "UMAP  (precomputed dist)"))
            elif m == "tsne":
                panels.append((dr.tsne_2d(D, precomputed=True), "t-SNE  (precomputed dist)"))
            elif m == "pcoa":
                c, (p1, p2) = dr.pcoa_2d(D); panels.append((c, f"PCoA  ({p1:.1f}%, {p2:.1f}%)"))
            else:
                raise SystemExit(f"method {m} not valid for --pairs")
    else:
        raise SystemExit("provide --features or --pairs")

    landscape_map.render_panels(panels, y, args.title, args.out,
                                footnote=args.footnote or None)
    print(f"saved {args.out}")


def _labels_for_features(df, ids, args):
    if args.label_col:
        return df[args.label_col].to_numpy()
    if args.labels_parallel:
        lab = pd.read_csv(args.labels_parallel)
        if len(lab) != len(df):
            raise SystemExit(f"--labels-parallel row count {len(lab)} != features {len(df)}")
        return lab[args.class_col].to_numpy()
    if args.labels_join:
        return _labels_join(ids, args.labels_join, args)
    raise SystemExit("need --label-col / --labels-parallel / --labels-join")


def _labels_join(ids, labels_csv, args):
    lab = pd.read_csv(labels_csv)
    id2c = lab.drop_duplicates(args.id_col).set_index(args.id_col)[args.class_col].to_dict()
    return np.array([id2c.get(str(i), -1) for i in ids])


if __name__ == "__main__":
    main()
