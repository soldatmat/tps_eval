from __future__ import annotations

import argparse

from domain_composition import extract_domain_composition_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract per-design TPS structural-domain composition (how many and "
        "which domain types EnzymeExplorer detects) from a directory of structures, "
        "writing a CSV keyed by ID for use as a filtration criterion alongside the other "
        "tps_eval metrics. EVERY input design gets one row, including designs with zero "
        "detected domains (n_domains=0). Requires the EnzymeExplorer env (PyMOL+foldseek; "
        "CPU-only)."
    )
    parser.add_argument(
        "structs_dir",
        help="Directory of generated structures; EE domain detection consumes the .pdb "
        "files (ID = filename stem).",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_domain_composition.csv next to the dir).",
    )
    parser.add_argument(
        "--detections_json",
        default=None,
        help="Path to an existing EE domain-detection JSON sidecar to PARSE instead of "
        "re-running detection (cheap re-use). If given but missing, detection runs and "
        "writes here.",
    )
    parser.add_argument("--n_jobs", type=int, default=10, help="Parallel jobs for detection.")
    parser.add_argument("--n_iters", type=int, default=3, help="EE detection iterations.")
    args = parser.parse_args()

    extract_domain_composition_dir(
        args.structs_dir,
        save_path=args.save_path,
        detections_json=args.detections_json,
        n_jobs=args.n_jobs,
        n_iters=args.n_iters,
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
