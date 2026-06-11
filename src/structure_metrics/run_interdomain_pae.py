from __future__ import annotations

import argparse

from interdomain_pae import extract_interdomain_pae_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inter-domain PAE: the AlphaFold/ESMFold Predicted-Aligned-Error "
        "confidence in the RELATIVE orientation of a design's TPS domains — a "
        "multi-domain packing failure mode per-residue pLDDT cannot see. Per design, "
        "gets per-domain residue ranges from the EnzymeExplorer domain detector "
        "(same detect_domains as domain_composition, keeping the residue spans), "
        "loads the <ID>_pae.npz saved at fold time, and averages the off-diagonal "
        "inter-domain PAE blocks (both directions). Writes a CSV keyed by ID. "
        "Single-/zero-domain designs and designs without a PAE npz are N/A (NaN). "
        "Requires the EnzymeExplorer env (PyMOL+foldseek; CPU-only)."
    )
    parser.add_argument(
        "structs_dir",
        help="Directory of generated structures; EE domain detection consumes the "
        ".pdb files (ID = filename stem).",
    )
    parser.add_argument(
        "--pae_dir",
        required=True,
        help="Directory of <ID>_pae.npz matrices saved at fold time (esmfold.py "
        "--pae_dir / the AF3 PAE extractor). ID = structure filename stem.",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_interdomain_pae.csv next to the dir).",
    )
    parser.add_argument(
        "--detections_json",
        default=None,
        help="Path to an existing EE domain-detection JSON sidecar to PARSE instead "
        "of re-running detection (cheap re-use; e.g. one written by domain_composition). "
        "If given but missing, detection runs and writes here.",
    )
    parser.add_argument(
        "--per_pair",
        action="store_true",
        help="Also emit a pae_<A>_<B> column per inter-domain pair (per design).",
    )
    parser.add_argument("--n_jobs", type=int, default=10, help="Parallel jobs for detection.")
    parser.add_argument("--n_iters", type=int, default=3, help="EE detection iterations.")
    args = parser.parse_args()

    extract_interdomain_pae_dir(
        args.structs_dir,
        args.pae_dir,
        save_path=args.save_path,
        detections_json=args.detections_json,
        per_pair=args.per_pair,
        n_jobs=args.n_jobs,
        n_iters=args.n_iters,
    )


if __name__ == "__main__":
    main()
