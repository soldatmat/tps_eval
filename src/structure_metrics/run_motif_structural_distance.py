from __future__ import annotations

import argparse

from motif_structural_distance import motif_structural_distance_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="3D (Angstrom) distance between the two class-I TPS metal-binding "
        "motifs (DDXXD-family and NSE/DTE) for every structure in a directory. "
        "Fold-agnostic: works for AlphaFold and ESMFold .pdb/.cif alike. Derives the "
        "sequence from each structure, localizes both motifs, and measures the "
        "centroid (and min CA-CA) distance between their metal-coordinating residues. "
        "Writes a CSV keyed by ID (filename stem); distances are NaN when a motif is "
        "absent."
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
        help="Output CSV path (default: <structs_dir>_motif_structural_distance.csv "
        "next to the directory).",
    )
    args = parser.parse_args()
    motif_structural_distance_dir(args.structs_dir, save_path=args.save_path)


if __name__ == "__main__":
    main()
