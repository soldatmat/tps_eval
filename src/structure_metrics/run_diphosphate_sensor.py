from __future__ import annotations

import argparse

from diphosphate_sensor import (
    DEFAULT_CUTOFF,
    DEFAULT_RY_DIST,
    diphosphate_sensor_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Diphosphate-sensor active-site metric for class-I TPS designs: "
        "does the design supply the basic residues (Arg/Lys) and the conserved RY pair "
        "that anchor and help ionize the substrate's diphosphate at the metal site? "
        "Anchors on the carboxylate-cage metal point (DDXXD-only fallback), counts the "
        "Arg/Lys whose terminal N atoms are near it AND point toward it, and detects "
        "sensor-Arg / Tyr (RY) pairs. Fold-agnostic (AlphaFold/ESMFold .pdb/.cif), "
        "apo-robust (no metals/substrate needed). Writes a CSV keyed by ID (filename "
        "stem); metal_point_found=False with zero counts when the site can't be located."
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
        help="Output CSV path (default: <structs_dir>_diphosphate_sensor.csv next to "
        "the directory).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF,
        help=f"Distance (A) from the metal point within which a basic residue's "
        f"terminal N must lie to count (default {DEFAULT_CUTOFF}).",
    )
    parser.add_argument(
        "--ry_dist",
        type=float,
        default=DEFAULT_RY_DIST,
        help=f"Distance (A) for the spatial RY-pair criterion: a Tyr OH within this of "
        f"a sensor Arg's guanidinium centroid (default {DEFAULT_RY_DIST}).",
    )
    args = parser.parse_args()
    diphosphate_sensor_dir(
        args.structs_dir,
        save_path=args.save_path,
        cutoff=args.cutoff,
        ry_dist=args.ry_dist,
    )


if __name__ == "__main__":
    main()
