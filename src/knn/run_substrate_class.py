"""argv entry for the substrate-class combiner.

Fuses the k-NN substrate vote (the three similarity-space top-k CSVs + the substrate
label file + the substrate calibration JSON) with the pocket-volume band and the
EnzymeExplorer sequence-only substrate signal into one <input>_substrate_class.csv.

At least one of the three --*_topk CSVs is required (the k-NN call is the primary signal);
--pocket_csv and --ee_csv are optional corroborating signals.

Example:
    python run_substrate_class.py \
        --sequence_topk gen_max_sequence_identity_topk.csv \
        --embedding_topk gen_embedding_esm1b_min_embedding_distance_topk.csv \
        --structural_topk structs_structural_identity_topk.csv \
        --label_file substrate_labels.csv \
        --calibration knn_calibration_substrate.json \
        --pocket_csv structs_pocket_descriptors.csv \
        --ee_csv gen_enzyme_explorer_sequence_only.csv \
        --out gen_substrate_class.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from knn_label_transfer import load_calibration  # noqa: E402
from substrate_class import combine_substrate_class  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sequence_topk", default=None,
                    help="<input>_max_sequence_identity_topk.csv (score = identity %).")
    ap.add_argument("--embedding_topk", default=None,
                    help="<input>_min_embedding_distance_topk.csv (score = L2 distance).")
    ap.add_argument("--structural_topk", default=None,
                    help="<structs_dir>_structural_identity_topk.csv (score = TM-score).")
    ap.add_argument("--label_file", required=True,
                    help="SUBSTRATE reference_id,label CSV (substrate_labels.csv).")
    ap.add_argument("--calibration", required=True,
                    help="Substrate calibration JSON (knn_calibration_substrate.json).")
    ap.add_argument("--pocket_csv", default=None,
                    help="<structs_dir>_pocket_descriptors.csv (catalytic_pocket_volume).")
    ap.add_argument("--ee_csv", default=None,
                    help="<input>_enzyme_explorer_sequence_only.csv (per-substrate scores).")
    ap.add_argument("--top_k", type=int, default=None,
                    help="Cap neighbours per query (default: all present).")
    ap.add_argument("--out", required=True, help="Output predictions CSV (keyed by ID).")
    args = ap.parse_args()

    spaces = {}
    if args.sequence_topk:
        spaces["sequence"] = args.sequence_topk
    if args.embedding_topk:
        spaces["embedding"] = args.embedding_topk
    if args.structural_topk:
        spaces["structural"] = args.structural_topk
    if not spaces:
        raise SystemExit(
            "Provide at least one of --sequence_topk / --embedding_topk / --structural_topk."
        )

    calibration = load_calibration(args.calibration)
    df = combine_substrate_class(
        spaces,
        args.label_file,
        calibration,
        pocket_csv=args.pocket_csv,
        ee_csv=args.ee_csv,
        top_k=args.top_k,
    )
    df.to_csv(args.out, index=False)

    n_unknown = int((df["predicted_substrate"] == "unknown").sum())
    n_pocket = int(df["substrate_agreement"].astype(str).eq("True").sum())
    n_ee = int(df["ee_agreement"].astype(str).eq("True").sum())
    print(f"Wrote {len(df)} substrate-class predictions to {args.out} "
          f"({n_unknown} unknown).")
    print(f"  pocket-volume band agrees with k-NN: {n_pocket}/{len(df)}")
    print(f"  EnzymeExplorer argmax agrees with k-NN: {n_ee}/{len(df)}")
    print("  predicted_substrate distribution:",
          df["predicted_substrate"].value_counts().to_dict())


if __name__ == "__main__":
    main()
