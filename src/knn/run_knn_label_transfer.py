"""argv entry for the k-NN coarse-label transfer tool.

Two subcommands:

  calibrate   Leave-one-out calibration on MARTS-DB self top-k CSVs -> JSON artifact.
  predict     Transfer labels to designs from their top-k CSVs + a calibration JSON.

Both are label-agnostic: the labeling is whatever the --label_file maps. The three
spaces are passed by their top-k CSV paths; a space may be omitted (the design/query
simply abstains there).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from knn_label_transfer import (  # noqa: E402
    calibrate,
    load_calibration,
    save_calibration,
    transfer_labels,
)


def _collect_spaces(args) -> dict:
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
    return spaces


def _add_space_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--sequence_topk", default=None,
                   help="<input>_max_sequence_identity_topk.csv (score = identity %).")
    p.add_argument("--embedding_topk", default=None,
                   help="<input>_min_embedding_distance_topk.csv (score = L2 distance).")
    p.add_argument("--structural_topk", default=None,
                   help="<structs_dir>_structural_identity_topk.csv (score = TM-score).")
    p.add_argument("--label_file", required=True,
                   help="CSV mapping reference_id,label (the labeling is the input).")
    p.add_argument("--top_k", type=int, default=None,
                   help="Cap neighbours per query (default: all present in the CSVs).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("calibrate", help="LOO calibration on MARTS-DB self top-k CSVs.")
    _add_space_args(pc)
    pc.add_argument("--out", required=True, help="Output calibration JSON path.")
    pc.add_argument("--labeling", default="labeling",
                    help="Name recorded in the artifact (e.g. first_cyclization).")
    pc.add_argument("--target_accuracy", type=float, default=0.5,
                    help="Accuracy floor used to pick tau per space/ensemble (default 0.5).")

    pp = sub.add_parser("predict", help="Transfer labels to designs.")
    _add_space_args(pp)
    pp.add_argument("--calibration", required=True, help="Calibration JSON from `calibrate`.")
    pp.add_argument("--out", required=True, help="Output predictions CSV (keyed by ID).")

    args = parser.parse_args()
    spaces = _collect_spaces(args)

    if args.cmd == "calibrate":
        cal = calibrate(
            spaces,
            args.label_file,
            top_k=args.top_k,
            target_accuracy=args.target_accuracy,
            labeling=args.labeling,
        )
        save_calibration(cal, args.out)
        print(f"Wrote calibration to {args.out}")
        print(f"  labeling={cal['labeling']}  n_classes={cal['n_classes']}")
        for space, s in cal["spaces"].items():
            print(
                f"  [{space:10s}] tau={s['tau']:.3f}  acc={s['accuracy']:.3f}  "
                f"predicted={s['n_predicted']}/{s['n_queries']}  abstained={s['n_abstained']}"
            )
        e = cal["ensemble"]
        print(
            f"  [ensemble  ] acc={e['accuracy']:.3f}  "
            f"predicted={e['n_predicted']}/{e['n_queries']}  abstained={e['n_abstained']}"
        )
    elif args.cmd == "predict":
        cal = load_calibration(args.calibration)
        df = transfer_labels(spaces, args.label_file, cal, top_k=args.top_k)
        df.to_csv(args.out, index=False)
        n_abstain = int((df["predicted_label"] == "unknown").sum())
        print(f"Wrote {len(df)} predictions to {args.out} ({n_abstain} abstained).")


if __name__ == "__main__":
    main()
