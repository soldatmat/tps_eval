from __future__ import annotations

import argparse

from radius_of_gyration import radius_of_gyration_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Radius-of-gyration / compactness structural metrics for TPS "
        "designs. Over the Cα atoms of each structure (all protein chains of the "
        "first model), emits RAW geometric numbers only: the unweighted "
        "radius_of_gyration (Å), asphericity and acylindricity from the gyration "
        "tensor, the three principal radii of gyration, and the Cα count. No "
        "compactness ratio or expected-Rg band is computed (that comparison is done "
        "downstream). Fold-agnostic: works for AlphaFold and ESMFold .pdb/.cif. "
        "Writes a CSV keyed by ID; unparsable structures yield a NaN row."
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
        help="Output CSV path (default: <structs_dir>_radius_of_gyration.csv "
        "next to the directory).",
    )
    args = parser.parse_args()
    radius_of_gyration_dir(args.structs_dir, save_path=args.save_path)


if __name__ == "__main__":
    main()
