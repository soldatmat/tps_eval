from __future__ import annotations

import argparse

from substrate_positioning import (
    DEFAULT_COORD_CUTOFF,
    DEFAULT_ION_RESNAMES,
    DEFAULT_MIN_SUBSTRATE_CARBONS,
    DEFAULT_SITE_RADIUS,
    substrate_positioning_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Substrate-positioning check for AF3 holo-folded class-I TPS designs: is "
        "the co-folded prenyl-diphosphate substrate (--af3_cofold mg_<sub> | mg_ee) poised in "
        "the carboxylate cage -- its diphosphate at the DDXXD/NSE metal cage and its reactive "
        "carbon (C1) held near the catalytic machinery? Auto-detects the substrate ligand by "
        "composition (>=1 P + carbons), so it works for forced and per-design co-folds alike; "
        "the bare POP pyrophosphate of mg_ppi and the Mg/Mn ions are excluded. Fold-agnostic "
        "(AlphaFold/ESMFold .pdb/.cif); apo / Mg-only / mg_ppi structures carry no substrate -> "
        "a graceful not-applicable row (substrate_present=False). Writes a CSV keyed by ID."
    )
    parser.add_argument(
        "structs_dir",
        help="Either an AlphaFold3 af_output directory (per-job subfolders with "
        "<job>/<job>_model.cif; ID = job name) OR a flat directory of .pdb/.cif structures "
        "(ID = filename stem).",
    )
    parser.add_argument(
        "--save_path", default=None,
        help="Output CSV path (default: <structs_dir>_substrate_positioning.csv next to the dir).",
    )
    parser.add_argument(
        "--site_radius", type=float, default=DEFAULT_SITE_RADIUS,
        help=f"Reported diphosphate-centroid -> cage-centroid distance band (informational; "
        f"default {DEFAULT_SITE_RADIUS}).",
    )
    parser.add_argument(
        "--coord_cutoff", type=float, default=DEFAULT_COORD_CUTOFF,
        help=f"Diphosphate-atom -> cage-oxygen distance (A) under which substrate_in_site is "
        f"True (the robust in-cage test; default {DEFAULT_COORD_CUTOFF}).",
    )
    parser.add_argument(
        "--ion_resnames", nargs="+", default=list(DEFAULT_ION_RESNAMES),
        help=f"Ion HETATM residue names to exclude from substrate detection and to measure the "
        f"diphosphate->ion distance against (default {' '.join(DEFAULT_ION_RESNAMES)}).",
    )
    parser.add_argument(
        "--min_substrate_carbons", type=int, default=DEFAULT_MIN_SUBSTRATE_CARBONS,
        help=f"Min carbons (with >=1 P) for a HETATM residue to count as a prenyl-PP substrate "
        f"(default {DEFAULT_MIN_SUBSTRATE_CARBONS}).",
    )
    parser.add_argument(
        "--substrate_resname", default=None,
        help="Force a specific ligand residue name as the substrate (default: auto-detect by "
        "composition).",
    )
    args = parser.parse_args()
    substrate_positioning_dir(
        args.structs_dir,
        save_path=args.save_path,
        site_radius=args.site_radius,
        coord_cutoff=args.coord_cutoff,
        ion_resnames=tuple(args.ion_resnames),
        min_substrate_carbons=args.min_substrate_carbons,
        substrate_resname=args.substrate_resname,
    )


if __name__ == "__main__":
    main()
