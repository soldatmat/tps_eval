from __future__ import annotations

import argparse

from ion_site_check import (
    DEFAULT_COORD_CUTOFF,
    DEFAULT_DIPHOSPHATE_RESNAMES,
    DEFAULT_ION_RESNAMES,
    DEFAULT_MIN_COORD_CONTACTS,
    DEFAULT_SITE_RADIUS,
    ion_site_check_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ion-placement (catalytic-site) check for AF3 holo-folded class-I "
        "TPS designs: does the trinuclear Mg2+/Mn2+ cluster AF3 co-folded "
        "(--af3_cofold mg|mg_ppi) actually land in the carboxylate cage? Unlike the "
        "apo-robust active-site tools, this READS the ion HETATMs and compares them to "
        "the expected apo cage point (active_site_geometry.metal_point): distance to "
        "the cage centroid, ions inside the site sphere, and per-ion Mg-O coordination "
        "contacts. Fold-agnostic (AlphaFold/ESMFold .pdb/.cif); apo / ESMFold "
        "structures carry no ions -> a graceful not-applicable row (n_ions_modelled=0). "
        "Writes a CSV keyed by ID (filename stem)."
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
        help="Output CSV path (default: <structs_dir>_ion_site_check.csv next to the "
        "directory).",
    )
    parser.add_argument(
        "--site_radius",
        type=float,
        default=DEFAULT_SITE_RADIUS,
        help=f"Distance (A) from the cage centroid within which an ion counts as "
        f"in-site (default {DEFAULT_SITE_RADIUS}).",
    )
    parser.add_argument(
        "--coord_cutoff",
        type=float,
        default=DEFAULT_COORD_CUTOFF,
        help=f"Mg-O coordination distance cutoff (A); real bonds ~2.0-2.5 A "
        f"(default {DEFAULT_COORD_CUTOFF}).",
    )
    parser.add_argument(
        "--min_coord_contacts",
        type=int,
        default=DEFAULT_MIN_COORD_CONTACTS,
        help=f"Min coordinating-oxygen contacts for an ion to count as coordinated / "
        f"well-placed (default {DEFAULT_MIN_COORD_CONTACTS}).",
    )
    parser.add_argument(
        "--ion_resnames",
        nargs="+",
        default=list(DEFAULT_ION_RESNAMES),
        help=f"Ion HETATM residue names to read (default {' '.join(DEFAULT_ION_RESNAMES)}).",
    )
    parser.add_argument(
        "--diphosphate_resnames",
        nargs="+",
        default=list(DEFAULT_DIPHOSPHATE_RESNAMES),
        help=f"Diphosphate HETATM residue names to read for the mg_ppi case "
        f"(default {' '.join(DEFAULT_DIPHOSPHATE_RESNAMES)}).",
    )
    args = parser.parse_args()
    ion_site_check_dir(
        args.structs_dir,
        save_path=args.save_path,
        site_radius=args.site_radius,
        coord_cutoff=args.coord_cutoff,
        min_coord_contacts=args.min_coord_contacts,
        ion_resnames=tuple(args.ion_resnames),
        diphosphate_resnames=tuple(args.diphosphate_resnames),
    )


if __name__ == "__main__":
    main()
