from __future__ import annotations

"""argv entry for the specificity-determining-residue (SDR) divergence tool.

Flags each design that is GLOBALLY close to a known-product TPS but DIVERGES at the
active-site residues that determine product specificity (the TEAS/HPS single-switch
failure mode global-similarity transfer misses). Consumes the committed --top_k
neighbour CSVs to pick the nearest known-TPS neighbour, then structurally aligns and
compares the SDR panel. See sdr_divergence.py for the method and the panel format.
"""

import argparse
import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from sdr_divergence import (  # noqa: E402
    DEFAULT_MAP_TOLERANCE,
    DEFAULT_PANEL_CUTOFF,
    DEFAULT_TAU_HIGH,
    DEFAULT_TAU_LOW,
    sdr_divergence_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "structs_dir",
        help="Directory of design structures (AF3 af_output or flat .pdb/.cif; ID = stem).",
    )
    parser.add_argument(
        "known_structs_dir",
        help="Directory of known-TPS reference structures the neighbour ids resolve to.",
    )
    parser.add_argument(
        "--structural_topk", default=None,
        help="<structs_dir>_structural_identity_topk.csv (rank-1 = nearest by TM-score; preferred).",
    )
    parser.add_argument(
        "--sequence_topk", default=None,
        help="<input>_max_sequence_identity_topk.csv (rank-1 = nearest by identity %; fallback).",
    )
    parser.add_argument(
        "--sdr_panel", default=None,
        help="Optional explicit SDR-panel CSV (reference_id,resnum[,expected_residue]) "
        "anchored to named references. Default: structure-derived active-site panel.",
    )
    parser.add_argument(
        "--panel_cutoff", type=float, default=DEFAULT_PANEL_CUTOFF,
        help=f"Structure-derived panel: heavy-atom cutoff (A) around the metal point "
        f"(default {DEFAULT_PANEL_CUTOFF}).",
    )
    parser.add_argument(
        "--map_tolerance", type=float, default=DEFAULT_MAP_TOLERANCE,
        help=f"Max Ca-Ca distance (A) to accept a design residue as a neighbour SDR "
        f"residue's counterpart after superposition (default {DEFAULT_MAP_TOLERANCE}).",
    )
    parser.add_argument(
        "--tau_high", type=float, default=DEFAULT_TAU_HIGH,
        help=f"Global-similarity floor in [0,1] for 'looks like this known TPS' "
        f"(default {DEFAULT_TAU_HIGH}).",
    )
    parser.add_argument(
        "--tau_low", type=float, default=DEFAULT_TAU_LOW,
        help=f"SDR-identity ceiling in [0,1] for 'differs where it matters' "
        f"(default {DEFAULT_TAU_LOW}).",
    )
    parser.add_argument(
        "--save_path", default=None,
        help="Output CSV path (default <structs_dir>_sdr_divergence.csv next to the dir).",
    )
    args = parser.parse_args()

    sdr_divergence_dir(
        args.structs_dir,
        args.known_structs_dir,
        structural_topk=args.structural_topk,
        sequence_topk=args.sequence_topk,
        panel_file=args.sdr_panel,
        panel_cutoff=args.panel_cutoff,
        map_tolerance=args.map_tolerance,
        tau_high=args.tau_high,
        tau_low=args.tau_low,
        save_path=args.save_path,
    )


if __name__ == "__main__":
    main()
