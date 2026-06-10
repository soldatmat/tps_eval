from __future__ import annotations

import argparse

from active_site_geometry import active_site_geometry_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Active-site-specific structural metrics for class-I TPS designs: "
        "the side-chain carboxylate-cage geometry of the metal-binding site "
        "(convergence radius of the coordinating oxygens, their count, and the "
        "metal-point void clearance at their centroid), plus an optional "
        "catalytic-constellation RMSD against reference templates. Fold-agnostic: "
        "works for AlphaFold and ESMFold .pdb/.cif alike. Derives the sequence from "
        "each structure, localizes the DDXXD-family and NSE/DTE motifs, and measures "
        "the catalytic site geometry. Writes a CSV keyed by ID (filename stem); "
        "geometry is NaN when a motif is absent."
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
        help="Output CSV path (default: <structs_dir>_active_site_geometry.csv "
        "next to the directory).",
    )
    parser.add_argument(
        "--templates",
        default=None,
        help="Comma-separated reference IDs (filename stems) in structs_dir to build "
        "the catalytic-constellation template from (e.g. '1ps1,5eat'). When given, "
        "the catalytic_constellation_rmsd / best_template columns are populated; "
        "otherwise they are NaN.",
    )
    args = parser.parse_args()
    template_ids = [t.strip() for t in args.templates.split(",") if t.strip()] if args.templates else None
    active_site_geometry_dir(args.structs_dir, save_path=args.save_path, template_ids=template_ids)


if __name__ == "__main__":
    main()
