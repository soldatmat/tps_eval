"""argv entry for the cyclization-geometry structure tool. See cyclization_geometry.py."""
import argparse

from cyclization_geometry import (
    DEFAULT_AROMATIC_CUTOFF,
    DEFAULT_FARCHAIN_BONDS,
    cyclization_geometry_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cyclization-relevant holo geometry for AF3-cofolded class-I TPS designs: "
        "substrate fold (rgyr, C1->distal fold-back, end-to-end) + aromatic cation-pi track. "
        "Necessary-not-sufficient cyclization signals; reference-independent (no apo metal_point). "
        "Apo / no-substrate folds -> graceful not-applicable rows. Writes a CSV keyed by ID."
    )
    parser.add_argument(
        "structs_dir",
        help="An AlphaFold3 af_output directory (per-job <job>/<job>_model.cif; ID = job name) "
        "OR a flat directory of .pdb/.cif structures (ID = filename stem).",
    )
    parser.add_argument(
        "--save_path", default=None,
        help="Output CSV path (default: <structs_dir>_cyclization_geometry.csv next to the dir).",
    )
    parser.add_argument(
        "--aromatic_cutoff", type=float, default=DEFAULT_AROMATIC_CUTOFF,
        help=f"Substrate-carbon -> aromatic-ring-centroid distance (A) for a cation-pi contact "
        f"(default {DEFAULT_AROMATIC_CUTOFF}).",
    )
    parser.add_argument(
        "--farchain_bonds", type=int, default=DEFAULT_FARCHAIN_BONDS,
        help=f"A chain carbon this many C-C bonds (or more) from C1 counts as distal for the "
        f"fold-back test (default {DEFAULT_FARCHAIN_BONDS}).",
    )
    parser.add_argument(
        "--ion_resnames", nargs="+", default=["MG", "MN"],
        help="Ion HETATM residue names to exclude from substrate detection (default MG MN).",
    )
    parser.add_argument(
        "--min_substrate_carbons", type=int, default=5,
        help="Min carbons (with >=1 P) for a HETATM residue to count as a prenyl-PP substrate "
        "(default 5).",
    )
    parser.add_argument(
        "--substrate_resname", default=None,
        help="Force a specific ligand residue name as the substrate (default: auto-detect).",
    )
    args = parser.parse_args()
    cyclization_geometry_dir(
        args.structs_dir,
        save_path=args.save_path,
        aromatic_cutoff=args.aromatic_cutoff,
        farchain_bonds=args.farchain_bonds,
        ion_resnames=tuple(args.ion_resnames),
        min_substrate_carbons=args.min_substrate_carbons,
        substrate_resname=args.substrate_resname,
    )


if __name__ == "__main__":
    main()
