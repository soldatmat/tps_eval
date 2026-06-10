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
    main()
