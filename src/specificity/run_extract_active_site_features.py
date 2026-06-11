from __future__ import annotations

"""argv entry for active-site / cation-specific-residue feature extraction.

Computes the per-protein active-site feature vector (see
``extract_active_site_features.py``) for the MARTS-DB ESMFold structures and
merges in the first-cyclization product-class label so the downstream UMAP map
can colour points by class.

Usage:
    python run_extract_active_site_features.py \
        --structs_dir <dir of marts_E*.pdb> \
        --marts_csv  <TPS_first_cyclization.csv> \
        --out        <active_site_features.csv> \
        [--radius 12.0]

The output CSV is keyed by ``id`` (== Enzyme_marts_ID == structure stem) with:
    id, product_class_id, product_class_marts_id, substrate_name,
    metal_point_found, n_shell_residues, n_residues, radius_A,
    <32 feature columns>

Only enzymes present in the MARTS-DB first-cyclization CSV are emitted (the
structure dir contains extra MARTS-DB sequences not in the first-cyclization
reference set; those are dropped). Enzymes mapping to multiple product classes
keep the FIRST (lowest class id) as the primary colour label, with the full
set in ``all_product_class_ids`` (semicolon-joined) for transparency.

OSC / class-12 outlier
----------------------
Product class 12 = oxidosqualene cyclases (substrate (S)-2,3-epoxysqualene): a
DIFFERENT fold, no Mg2+ cluster, no real DDXXD/NSE-DTE motif. CAVEAT: the shared
relaxed DDXXD regex ([DE][DE]..[DE]) is permissive enough to spuriously match an
N-terminal acidic stretch in OSC sequences, so ``metal_point_found`` does NOT by
itself flag class-12 (its feature row is anchored on a non-catalytic motif and is
therefore NOT meaningful). We exclude class-12 by LABEL by default
(``--exclude_class 12``); pass ``--exclude_class ""`` to keep it (flagged).
"""

import argparse
import os
import sys

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from extract_active_site_features import (  # noqa: E402
    DEFAULT_RADIUS,
    active_site_features_dir,
)


def _class_metadata(marts_csv: str) -> pd.DataFrame:
    """One row per Enzyme_marts_ID with its product-class label(s)."""
    df = pd.read_csv(marts_csv)
    needed = {"Enzyme_marts_ID", "First_cyclization_product_id"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{marts_csv} missing column(s): {sorted(missing)}")

    sub_col = "Substrate_name" if "Substrate_name" in df.columns else None
    prod_col = "First_cyclization_product" if "First_cyclization_product" in df.columns else None

    rows = []
    for eid, grp in df.groupby("Enzyme_marts_ID"):
        classes = sorted(int(c) for c in grp["First_cyclization_product_id"].dropna().unique())
        primary = classes[0] if classes else None
        # marts product id of the primary class
        prim_marts = ""
        if prod_col is not None and primary is not None:
            m = grp.loc[grp["First_cyclization_product_id"] == primary, prod_col]
            if len(m):
                prim_marts = str(m.iloc[0])
        substrate = ""
        if sub_col is not None and primary is not None:
            s = grp.loc[grp["First_cyclization_product_id"] == primary, sub_col]
            if len(s):
                substrate = str(s.iloc[0])
        rows.append(
            {
                "id": eid,
                "product_class_id": primary,
                "product_class_marts_id": prim_marts,
                "substrate_name": substrate,
                "n_product_classes": len(classes),
                "all_product_class_ids": ";".join(str(c) for c in classes),
            }
        )
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structs_dir", required=True, help="Dir of marts_E*.pdb ESMFold structures.")
    p.add_argument("--marts_csv", required=True, help="TPS_first_cyclization.csv path.")
    p.add_argument("--out", required=True, help="Output CSV path.")
    p.add_argument("--radius", type=float, default=DEFAULT_RADIUS, help="Active-site shell radius (A).")
    p.add_argument(
        "--exclude_class",
        default="12",
        help="Comma-separated product-class ids to drop by LABEL (default '12' = "
        "the OSC outlier, whose active-site frame is incommensurable with class-I "
        "TPS). Pass '' to keep all classes.",
    )
    p.add_argument(
        "--features_only_csv",
        default=None,
        help="Optional path to also dump the raw structure-keyed feature CSV "
        "(before the MARTS-DB class merge / filter).",
    )
    args = p.parse_args(argv)

    meta = _class_metadata(args.marts_csv)
    exclude = {int(x) for x in args.exclude_class.split(",") if x.strip() != ""}
    if exclude:
        n_before = len(meta)
        meta = meta[~meta["product_class_id"].isin(exclude)].reset_index(drop=True)
        print(f"Excluding product class(es) {sorted(exclude)} by label: "
              f"{n_before} -> {len(meta)} enzymes")
    enzyme_ids = set(meta["id"])
    print(f"MARTS-DB first-cyclization enzymes (after exclusion): {len(enzyme_ids)}")

    feat = active_site_features_dir(
        args.structs_dir,
        save_path=args.features_only_csv,
        radius=args.radius,
        id_filter=enzyme_ids,
    )

    merged = meta.merge(feat, on="id", how="inner")
    # column order: metadata first, then features
    meta_cols = [
        "id",
        "product_class_id",
        "product_class_marts_id",
        "substrate_name",
        "n_product_classes",
        "all_product_class_ids",
    ]
    feat_cols = [c for c in feat.columns if c != "id"]
    merged = merged[meta_cols + feat_cols].sort_values("id").reset_index(drop=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    merged.to_csv(args.out, index=False)

    n_with = int(merged["metal_point_found"].sum())
    n_without = len(merged) - n_with
    print(f"\nWrote {len(merged)} enzyme rows to {args.out}")
    print(f"  with metal point (class-I cage): {n_with}")
    print(f"  without (OSC class-12 / no cage): {n_without}")
    if "product_class_id" in merged:
        cov = (
            merged.assign(has=merged["metal_point_found"])
            .groupby("product_class_id")["has"]
            .agg(["size", "sum"])
        )
        print("\n  per-class coverage (class_id: n_total / n_with_features):")
        for cid, r in cov.iterrows():
            print(f"    {int(cid):2d}: {int(r['size']):4d} / {int(r['sum']):4d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
