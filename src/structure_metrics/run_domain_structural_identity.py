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
    parser.add_argument(
        "--self_mode",
        action="store_true",
        default=False,
        help="Searching a domain set against itself: drop hits to a reference domain "
        "originating from the query's own source structure before the best-hit "
        "reduction, so each design's best hit is its nearest OTHER known-TPS domain "
        "(leave-one-out) instead of the trivial self-match TM~1.0.",
    )
    args = parser.parse_args()

    extract_domain_structural_identity_dir(
        args.structs_dir,
        args.known_domain_structures_root,
        save_path=args.save_path,
        n_jobs=args.n_jobs,
        n_iters=args.n_iters,
        keep_detected_domains=args.keep_detected_domains,
        exclude_self=args.self_mode,
    )


if __name__ == "__main__":
    # EnzymeExplorer's in-process domain detection (PyMOL + a spawn-based
    # multiprocessing Pool) can leave helper threads / pool workers that
    # deadlock the interpreter's multiprocessing atexit teardown AFTER our
    # output CSV is already written and flushed -- observed as a ~20h hang of an
    # otherwise-finished run on a login node (under SLURM it would silently burn
    # walltime to TIMEOUT). main() writes and closes the CSV before returning,
    # so the result is safely on disk; bypass the hanging teardown with a hard
    # exit. A failure still surfaces a traceback + non-zero exit code first.
    import os
    import sys

    exit_code = 0
    try:
        main()
    except SystemExit as exc:
        exit_code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
    except BaseException:
        import traceback

        traceback.print_exc()
        exit_code = 1
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(exit_code)
