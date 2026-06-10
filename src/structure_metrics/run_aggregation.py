# -*- coding: utf-8 -*-
"""argv entry for the Aggrescan3D (A3D) aggregation metric. Runs in the
PYTHON 2.7 ``aggrescan3d`` conda env -- keep this file Py2-compatible."""

import argparse

from aggregation import extract_aggregation_dir


def main():
    parser = argparse.ArgumentParser(
        description="Run Aggrescan3D (STATIC mode) on every structure (.pdb/.cif) "
        "in a directory and write a CSV keyed by ID with per-structure "
        "aggregation-propensity scalars (a3d_avg/total/max/min/total_pos_score), "
        "for use as a filtration criterion alongside the other tps_eval metrics. "
        "A structure-based aggregation/expressibility signal orthogonal to the "
        "sequence-based SoluProt. Only static SASA-based scoring runs; the slow "
        "dynamic CABS-flex mode is never triggered."
    )
    parser.add_argument(
        "structs_dir",
        help="Either an AlphaFold3 af_output directory (per-job subfolders with "
        "<job>/<job>_model.cif; ID = job name) OR a flat directory of .pdb/.cif "
        "structures (ID = filename stem). .cif inputs are converted to PDB with "
        "Biopython before scoring.",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_aggregation.csv next to the directory).",
    )
    parser.add_argument(
        "--save_residue_scores",
        action="store_true",
        help="Also dump each structure's per-residue A3D score array to a side "
        "directory (one <ID>.csv per structure) for hotspot visualization. "
        "Default off.",
    )
    parser.add_argument(
        "--residue_scores_dir",
        default=None,
        help="Directory for the per-residue score files when --save_residue_scores "
        "is set (default: <structs_dir>_aggregation_residue_scores).",
    )
    args = parser.parse_args()
    extract_aggregation_dir(
        args.structs_dir,
        save_path=args.save_path,
        save_residue_scores=args.save_residue_scores,
        residue_scores_dir=args.residue_scores_dir,
    )


if __name__ == "__main__":
    main()
