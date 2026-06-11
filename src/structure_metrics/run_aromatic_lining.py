from __future__ import annotations

import argparse

from aromatic_lining import (
    DEFAULT_CATION_PI_MAX,
    DEFAULT_CATION_PI_MIN,
    DEFAULT_CUTOFF,
    DEFAULT_FACE_ANGLE_DEG,
    aromatic_lining_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aromatic / cation-pi pocket-lining metric for class-I TPS designs: "
        "counts and ring-orientation of the Trp/Tyr/Phe residues lining the catalytic "
        "pocket (a proxy for carbocation-stabilization / cyclization capability). "
        "Fold-agnostic (AlphaFold or ESMFold .pdb/.cif). Localizes the DDXXD-family + "
        "NSE/DTE motifs to place the carboxylate-cage metal point, selects pocket "
        "residues within a distance shell of it, and counts pocket aromatics + the "
        "subset whose ring face points at the cavity interior within cation-pi range. "
        "Writes a CSV keyed by ID (filename stem); counts are NaN and metal_point_found "
        "False when a motif / coordinating oxygen is absent."
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
        help="Output CSV path (default: <structs_dir>_aromatic_lining.csv next to the directory).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF,
        help=f"Pocket shell radius (A) around the metal point (default {DEFAULT_CUTOFF}).",
    )
    parser.add_argument(
        "--cation_pi_min",
        type=float,
        default=DEFAULT_CATION_PI_MIN,
        help=f"Min ring-centroid->locus distance for an inward-facing hit (default {DEFAULT_CATION_PI_MIN}).",
    )
    parser.add_argument(
        "--cation_pi_max",
        type=float,
        default=DEFAULT_CATION_PI_MAX,
        help=f"Max ring-centroid->locus distance for an inward-facing hit (default {DEFAULT_CATION_PI_MAX}).",
    )
    parser.add_argument(
        "--face_angle_deg",
        type=float,
        default=DEFAULT_FACE_ANGLE_DEG,
        help="Max angle (deg) between the ring normal and the centroid->locus vector "
        f"for a face-on (inward) hit (default {DEFAULT_FACE_ANGLE_DEG}).",
    )
    args = parser.parse_args()
    aromatic_lining_dir(
        args.structs_dir,
        save_path=args.save_path,
        cutoff=args.cutoff,
        cation_pi_min=args.cation_pi_min,
        cation_pi_max=args.cation_pi_max,
        face_angle_deg=args.face_angle_deg,
    )


if __name__ == "__main__":
    main()
