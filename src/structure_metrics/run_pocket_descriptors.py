from __future__ import annotations

import argparse

from pocket_descriptors import pocket_descriptors_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Active-site pocket descriptors for class-I TPS designs: geometric "
        "cavity descriptors from fpocket plus an ML ligandability cross-check from "
        "P2Rank, BOTH anchored on the active-site metal point (the centroid of the "
        "DDXXD + NSE/DTE coordinating side-chain oxygens). For each structure the "
        "catalytic pocket is the detected pocket enclosing / nearest that metal point; "
        "we emit its fpocket volume/hydrophobicity/enclosure/n-alpha-spheres (+ SASA & "
        "depth proxy) and its P2Rank ligandability score & rank. Fold-agnostic "
        "(AlphaFold/ESMFold .pdb/.cif). RAW numbers only; NaN columns when no pocket "
        "coincides with the metal point (a meaningful red flag). Writes a CSV keyed by "
        "ID (filename stem)."
    )
    parser.add_argument(
        "structs_dir",
        help="Either an AlphaFold3 af_output directory (per-job subfolders with "
        "<job>/<job>_model.cif; ID = job name) OR a flat directory of .pdb/.cif "
        "structures (ID = filename stem).",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_pocket_descriptors.csv next to "
        "the directory).",
    )
    parser.add_argument(
        "--fpocket_bin",
        default="fpocket",
        help="fpocket executable (default: 'fpocket' on PATH from the pocket env).",
    )
    parser.add_argument(
        "--p2rank_bin",
        default=None,
        help="P2Rank 'prank' executable (e.g. $P2RANK_PATH/prank). When omitted the "
        "P2Rank cross-check is skipped and its columns are NaN.",
    )
    args = parser.parse_args()
    pocket_descriptors_dir(
        args.structs_dir,
        save_path=args.save_path,
        fpocket_bin=args.fpocket_bin,
        p2rank_bin=args.p2rank_bin,
    )


if __name__ == "__main__":
    main()
