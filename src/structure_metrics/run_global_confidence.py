from __future__ import annotations

import argparse

from global_confidence import extract_global_confidence_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Global fold confidence: the AlphaFold/ESMFold pTM (predicted "
        "TM-score; 0-1, higher=better) WHOLE-FOLD confidence scalar, complementing "
        "the per-residue pLDDT. Reads the ptm (and iptm when present) saved at fold "
        "time into the <ID>_pae.npz files (esmfold.py / src/alphafold/"
        "extract_pae.py) from --pae_dir and writes a CSV keyed by ID. The iptm column "
        "appears only for multi-chain inputs; single-chain TPS get pTM only. NaN when "
        "the npz is missing or carries no ptm. Reads numpy npz only (TPS_EVAL_ENV; "
        "CPU-only, no structure parsing)."
    )
    parser.add_argument(
        "--pae_dir",
        required=True,
        help="Directory of <ID>_pae.npz files saved at fold time (esmfold.py "
        "--pae_dir / the AF3 PAE extractor). ID = npz filename stem.",
    )
    parser.add_argument(
        "--structs_dir",
        default=None,
        help="Optional structs dir used only to NAME the output CSV "
        "(<structs_dir>_global_confidence.csv), mirroring the structure branch. "
        "Defaults to naming off --pae_dir.",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir|pae_dir>_global_confidence.csv).",
    )
    args = parser.parse_args()

    extract_global_confidence_dir(
        args.pae_dir,
        structs_dir=args.structs_dir,
        save_path=args.save_path,
    )


if __name__ == "__main__":
    main()
