from __future__ import annotations

# TmProt melting-temperature (Tm) prediction — sequence-branch metric.
#
# TmProt (Loschmidt Laboratories) predicts a protein's melting temperature from
# sequence alone, using ESM-2 (650M) fine-tuned with a LoRA adapter. It is
# vendored at vendor/TmProt and its standalone CLI (`tmprot`) is installed into
# the TMPROT_ENV conda env via scripts/setup_tmprot.sh
# (`pip install -e vendor/TmProt/tmprot-1.0`).
#
# Implementation: shell out to the installed `tmprot` console command once for
# the whole FASTA (`tmprot -i <fasta> -o <tmp> -d ,`), then reshape its output
# CSV into the tps_eval convention — a CSV keyed by `ID` with a single RAW metric
# column `tm` (predicted Tm in degC). We deliberately drop TmProt's Rank and
# Thermostable columns: Thermostable is a threshold-dependent label and tps_eval
# tools emit RAW numbers only (bands/thresholds are applied downstream).
#
# TmProt skips sequences it cannot score (shorter than 20 AA, longer than 2000
# AA, or containing non-standard amino acids); those IDs are ABSENT from its
# output. We reindex to the full input FASTA ID set so every input gets a row,
# with NaN where TmProt produced no prediction. ID = FASTA record id (matches the
# other sequence-branch tools).

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers

COLUMNS = ["ID", "tm"]

# Native column in TmProt's output CSV holding the predicted melting temperature.
_TMPROT_TM_COLUMN = "Predicted Tm [°C]"
_TMPROT_ID_COLUMN = "ID"


def default_save_path(fasta_path: str, out_suffix: str = "tmprot") -> str:
    """<fasta_dir>/<fasta_stem>_<out_suffix>.csv — the sequence-branch naming rule."""
    fasta_path = os.fspath(fasta_path)
    directory = os.path.dirname(os.path.abspath(fasta_path))
    stem = os.path.splitext(os.path.basename(fasta_path))[0]
    return os.path.join(directory, f"{stem}_{out_suffix}.csv")


def run_tmprot_cli(
    fasta_path: str,
    out_dir: str,
    *,
    device: Optional[str] = None,
    tmprot_executable: str = "tmprot",
) -> str:
    """Run the installed `tmprot` CLI on a FASTA, returning the output CSV path.

    `device` optionally forces execution: "cpu" hides the GPU (the CLI otherwise
    auto-selects cuda:0 when available); "cuda"/None leave TmProt's default. The
    CLI writes <out_dir>/<fasta_stem>.csv, comma-delimited (we pass `-d ,`).
    Raises RuntimeError on a nonzero exit or a missing output file.
    """
    fasta_path = os.fspath(fasta_path)
    stem = os.path.splitext(os.path.basename(fasta_path))[0]
    cmd = [tmprot_executable, "-i", fasta_path, "-o", out_dir, "-d", ","]

    env = os.environ.copy()
    if device == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""

    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(
            f"tmprot failed on {fasta_path}:\n{proc.stdout}\n{proc.stderr}"
        )
    out_csv = os.path.join(out_dir, stem + ".csv")
    if not os.path.isfile(out_csv):
        raise RuntimeError(
            f"tmprot produced no output CSV {out_csv}\n{proc.stdout}\n{proc.stderr}"
        )
    return out_csv


def _read_tmprot_output(out_csv: str) -> dict:
    """Parse TmProt's output CSV into {ID: tm}. Empty (header-only) is allowed."""
    df = pd.read_csv(out_csv)
    if df.empty:
        return {}
    if _TMPROT_ID_COLUMN not in df.columns or _TMPROT_TM_COLUMN not in df.columns:
        raise RuntimeError(
            f"Unexpected tmprot output columns {list(df.columns)} in {out_csv}; "
            f"expected {_TMPROT_ID_COLUMN!r} and {_TMPROT_TM_COLUMN!r}."
        )
    return {
        str(row[_TMPROT_ID_COLUMN]): float(row[_TMPROT_TM_COLUMN])
        for _, row in df.iterrows()
    }


def score_fasta(
    fasta_path: str,
    *,
    save_path: Optional[str] = None,
    out_suffix: str = "tmprot",
    device: Optional[str] = None,
    tmprot_executable: str = "tmprot",
) -> pd.DataFrame:
    """Predict Tm for every sequence in a FASTA, writing a CSV keyed by ID.

    The output has one row per input FASTA record (reindexed to the full ID set,
    NaN where TmProt skipped the sequence), sorted by ID, with the RAW `tm` column.
    """
    fasta_path = os.fspath(fasta_path)
    records = load_fasta_sequences(fasta_path, load_identifiers=True)
    identifiers, _ = separate_identifiers(records)

    with tempfile.TemporaryDirectory(prefix="tmprot_") as tmp:
        out_csv = run_tmprot_cli(
            fasta_path, tmp, device=device, tmprot_executable=tmprot_executable
        )
        tm_by_id = _read_tmprot_output(out_csv)

    rows: List[dict] = [
        {"ID": identifier, "tm": tm_by_id.get(identifier, np.nan)}
        for identifier in identifiers
    ]
    df = pd.DataFrame(rows, columns=COLUMNS).sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = default_save_path(fasta_path, out_suffix=out_suffix)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
