from __future__ import annotations

import argparse

from plddt import CONFIDENT_THRESHOLD, extract_plddt_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract per-structure pLDDT summaries from AlphaFold structures "
        "(.pdb/.cif) in a directory, writing a CSV keyed by ID for use as a "
        "filtration criterion alongside the other tps_eval metrics. NOTE: reads "
        "pLDDT from the B-factor field, valid only for PREDICTED structures "
        "(AlphaFold/ESMFold) — experimental structures' B-factors are temperature "
        "factors, not confidence."
    )
    parser.add_argument(
        "structs_dir",
        help="Either an AlphaFold3 af_output directory (per-job subfolders with "
        "<job>/<job>_model.cif — pLDDT read from the authoritative top-ranked "
        "model; ID = job name) OR a flat directory of .pdb/.cif structures whose "
        "B-factor already holds pLDDT (e.g. AF2/ColabFold, or structs/ extracted "
        "by the patched cif_to_pdb; ID = filename stem). structs/*.pdb extracted "
        "by the OLD obabel converter are zeroed — re-extract or use af_output.",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_plddt.csv next to the directory).",
    )
    parser.add_argument(
        "--confident_threshold",
        type=float,
        default=CONFIDENT_THRESHOLD,
        help=f"pLDDT cutoff for frac_plddt_confident (default: {CONFIDENT_THRESHOLD}).",
    )
    args = parser.parse_args()
    extract_plddt_dir(
        args.structs_dir,
        save_path=args.save_path,
        confident_threshold=args.confident_threshold,
    )


if __name__ == "__main__":
    main()
