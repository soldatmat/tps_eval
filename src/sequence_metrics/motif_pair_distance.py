from __future__ import annotations

"""Sequence distance between the two class-I TPS metal-binding motifs.

For each sequence in a FASTA, locate the DDXXD-family motif and the NSE/DTE motif
(via the shared ``motif_localization`` helper) and report the residue separation
between them. Output is a CSV keyed by ``ID`` (one row per sequence), composable
with the other tps_eval sequence metrics. Distances are NaN when either motif is
absent.

Distance definitions (both signed positive when NSE/DTE follows DDXXD, the
canonical class-I order):

* ``motif_start_distance`` = ``start_NSE - start_DDXXD`` (1-based starts) — the
  separation between the two motifs' first residues. Negative if NSE/DTE precedes
  DDXXD.
* ``motif_gap`` = ``start_NSE - end_DDXXD`` (0-based) — the number of residues
  strictly between the end of DDXXD and the start of NSE/DTE (the inter-motif gap;
  negative if the motifs overlap).

Helper columns record each motif's matched substring and 1-based start so the
distance is auditable.
"""

import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers  # noqa: E402
from sequence_metrics.motif_localization import locate_ddxxd, locate_nse_dte  # noqa: E402

COLUMNS = [
    "ID",
    "motif_start_distance",
    "motif_gap",
    "ddxxd_motif",
    "ddxxd_start",
    "nse_dte_motif",
    "nse_dte_start",
]


def _get_save_path(fasta_path: str) -> str:
    extension = fasta_path.split(".")[-1]
    return fasta_path[: -len(extension) - 1] + "_motif_pair_distance.csv"


def _row(seq_id: str, sequence: str) -> dict:
    ddxxd = locate_ddxxd(sequence)
    nse = locate_nse_dte(sequence)
    if ddxxd is None or nse is None:
        start_distance: Optional[float] = np.nan
        gap: Optional[float] = np.nan
    else:
        start_distance = nse.start_1 - ddxxd.start_1   # 1-based start separation
        gap = nse.start - ddxxd.end                    # 0-based inter-motif gap
    return {
        "ID": seq_id,
        "motif_start_distance": start_distance,
        "motif_gap": gap,
        "ddxxd_motif": ddxxd.matched if ddxxd else "",
        "ddxxd_start": ddxxd.start_1 if ddxxd else np.nan,
        "nse_dte_motif": nse.matched if nse else "",
        "nse_dte_start": nse.start_1 if nse else np.nan,
    }


def motif_pair_distance(fasta_path: str, *, save: bool = True) -> pd.DataFrame:
    """Per-sequence residue distance between the DDXXD and NSE/DTE motifs."""
    identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    rows = [_row(str(i).split(" ", 1)[0], str(s)) for i, s in zip(identifiers, sequences)]
    df = pd.DataFrame(rows)[COLUMNS]

    if save:
        save_path = _get_save_path(fasta_path)
        df.to_csv(save_path, index=False)
        print(f"Wrote {len(df)} rows to {save_path}")
    return df
