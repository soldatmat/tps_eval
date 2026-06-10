from __future__ import annotations

import argparse

from proteinmpnn_score import score_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ProteinMPNN sequence likelihood (NLL) of each design's OWN "
        "sequence given its backbone, for every structure in a directory. "
        "Lower proteinmpnn_nll = the sequence is more compatible with the fold. "
        "Writes a CSV keyed by ID."
    )
    parser.add_argument(
        "structs_dir",
        help="AF3 af_output dir (<job>/<job>_model.cif) OR flat dir of .pdb/.cif; "
        "ID = job name / filename stem.",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_proteinmpnn_score.csv).",
    )
    parser.add_argument(
        "--model_name", default="v_48_020",
        help="ProteinMPNN model name (default: v_48_020).",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed (0 = random).")
    parser.add_argument(
        "--backbone_noise", type=float, default=0.0,
        help="Std of Gaussian noise added to backbone atoms when scoring (default 0).",
    )
    args = parser.parse_args()

    score_dir(
        args.structs_dir,
        save_path=args.save_path,
        model_name=args.model_name,
        seed=args.seed,
        backbone_noise=args.backbone_noise,
    )


if __name__ == "__main__":
    main()
