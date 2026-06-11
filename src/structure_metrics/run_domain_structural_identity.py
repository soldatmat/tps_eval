from __future__ import annotations

import argparse

from domain_structural_identity import extract_domain_structural_identity_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Domain-level structural identity: detect TPS domains in each "
        "generated structure (EnzymeExplorer) and foldseek-align them against the "
        "curated known-TPS reference DOMAINS, recording each design's best "
        "domain-level TM-score (+ lddt) to the nearest known-TPS domain, the matched "
        "reference domain + its type, and per-domain-type bests. The domain-level "
        "analog of structural_identity. Writes a CSV keyed by ID; designs with no "
        "detected domain get a NaN row. Requires an env with BOTH EnzymeExplorer "
        "(detect_domains) AND foldseek (the EE 'prod' env on Aurum has both)."
    )
    parser.add_argument(
        "structs_dir",
        help="Directory of generated structures; EE domain detection consumes the "
        ".pdb files (ID = filename stem).",
    )
    parser.add_argument(
        "known_domain_structures_root",
        help="Directory of known-TPS reference DOMAIN structures (EE's martsDB "
        "detected domains; flat dir of <stem>_<type>_<index>.pdb files).",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_domain_structural_identity.csv).",
    )
    parser.add_argument("--n_jobs", type=int, default=10, help="Parallel jobs for detection.")
    parser.add_argument("--n_iters", type=int, default=3, help="EE detection iterations.")
    parser.add_argument(
        "--keep_detected_domains",
        default=None,
        help="If given, write the per-design detected domain .pdb files into this "
        "directory and keep them (default: a temp dir, removed after).",
    )
    args = parser.parse_args()

    extract_domain_structural_identity_dir(
        args.structs_dir,
        args.known_domain_structures_root,
        save_path=args.save_path,
        n_jobs=args.n_jobs,
        n_iters=args.n_iters,
        keep_detected_domains=args.keep_detected_domains,
    )


if __name__ == "__main__":
    main()
