from __future__ import annotations

# CataPro enzyme-kinetics prediction — sequence-branch metric.
#
# For each sequence, predict the steady-state kinetic parameters of the reaction
# with a chosen substrate using CataPro (vendored at vendor/CataPro). CataPro is a
# ProtT5 (enzyme) + MolT5 (substrate) + MACCS-fingerprint model with 10 folds whose
# predictions are averaged; it reports three quantities on a log10 scale:
#   * pred_log10[kcat(s^-1)]            — turnover number
#   * pred_log10[Km(mM)]               — Michaelis constant
#   * pred_log10[kcat/Km(s^-1 mM^-1)]  — catalytic efficiency
# We exponentiate to absolute units and report catapro_kcat (s^-1), catapro_km (mM),
# catapro_kcat_km (s^-1 mM^-1). HIGHER kcat / kcat_km = more active on that substrate.
#
# SUBSTRATE: CataPro needs a substrate SMILES per sequence. TPS substrates are the
# prenyl-diphosphates in alphafold.cofold_substrates.SUBSTRATE_SMILES (the same map
# used for AF3 co-folding). We look up --target_substrate there. CataPro runs only
# where a SMILES is known; for a substrate code with no SMILES entry we emit a NaN
# row per the pipeline design ("runs where SMILES known, else NaN") rather than
# crashing. An explicit --smiles overrides the lookup.
#
# Implementation: build CataPro's input CSV (Enzyme_id, type=wild, sequence, smiles),
# shell out to vendor/CataPro/inference/predict.py (which MUST run with cwd set to
# that inference dir — it uses flat `from utils import *` imports), then reshape its
# output CSV to a tps_eval CSV keyed by ID. ID = FASTA header's first whitespace
# token; CataPro tags rows as "<ID>_wild", which we strip back to <ID>. Output is
# reindexed to the FULL FASTA ID set, so any sequence CataPro dropped gets a NaN row.

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers
from alphafold.cofold_substrates import SUBSTRATE_SMILES

REPO_ROOT = SRC_DIR.parent
CATAPRO_DIR = REPO_ROOT / "vendor" / "CataPro"
CATAPRO_INFERENCE_DIR = CATAPRO_DIR / "inference"
CATAPRO_PREDICT = CATAPRO_INFERENCE_DIR / "predict.py"
DEFAULT_MODEL_DPATH = CATAPRO_DIR / "models"

DEFAULT_OUT_SUFFIX = "catapro"

COLUMNS = ["ID", "catapro_kcat", "catapro_km", "catapro_kcat_km", "catapro_substrate"]

# Native (log10) column names written by vendor/CataPro/inference/predict.py.
NATIVE_ID = "fasta_id"
NATIVE_KCAT = "pred_log10[kcat(s^-1)]"
NATIVE_KM = "pred_log10[Km(mM)]"
NATIVE_KCAT_KM = "pred_log10[kcat/Km(s^-1mM^-1)]"

# CataPro tags every row's id as f"{Enzyme_id}_{type}"; we always submit type=wild.
CATAPRO_TYPE = "wild"
_ID_SUFFIX = f"_{CATAPRO_TYPE}"


def resolve_smiles(target_substrate: Optional[str], smiles: Optional[str] = None) -> Optional[str]:
    """Resolve the substrate SMILES to feed CataPro.

    An explicit ``smiles`` wins. Otherwise look ``target_substrate`` up in the
    shared SUBSTRATE_SMILES map (case-insensitive). Returns None when neither is
    available (caller then emits NaN rows).
    """
    if smiles:
        return smiles
    if not target_substrate:
        return None
    return SUBSTRATE_SMILES.get(target_substrate.upper())


def log10_to_absolute(log10_value: float) -> float:
    """Convert a CataPro log10 prediction to its absolute value (10 ** x)."""
    return float(10.0 ** log10_value)


def _default_save_path(fasta_path: str, out_suffix: str = DEFAULT_OUT_SUFFIX) -> str:
    """Sibling of the FASTA: <fasta_dir>/<fasta_stem>_<out_suffix>.csv."""
    directory = os.path.dirname(os.path.abspath(fasta_path))
    stem = os.path.splitext(os.path.basename(fasta_path))[0]
    return os.path.join(directory, f"{stem}_{out_suffix}.csv")


def _read_fasta_ids_and_sequences(fasta_path: str):
    """Return (ids, sequences) with ID = header's first whitespace token."""
    identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    ids = [str(identifier).split()[0] for identifier in identifiers]
    return ids, [str(sequence) for sequence in sequences]


def _nan_frame(ids: List[str], substrate: Optional[str]) -> pd.DataFrame:
    """All-NaN metric frame keyed by ID (used when no SMILES is available)."""
    frame = pd.DataFrame({"ID": ids})
    frame["catapro_kcat"] = np.nan
    frame["catapro_km"] = np.nan
    frame["catapro_kcat_km"] = np.nan
    frame["catapro_substrate"] = substrate
    return frame[COLUMNS].sort_values("ID").reset_index(drop=True)


def run_catapro_predict(
    input_csv: str,
    output_csv: str,
    *,
    model_dpath: str,
    device: str = "cuda",
    batch_size: int = 64,
    python_exe: str = sys.executable,
) -> None:
    """Shell out to CataPro's predict.py (cwd = its inference dir; flat imports).

    Raises RuntimeError on a nonzero exit or a missing output file.
    """
    cmd = [
        python_exe,
        str(CATAPRO_PREDICT),
        "-inp_fpath", os.path.abspath(input_csv),
        "-model_dpath", os.path.abspath(model_dpath),
        "-batch_size", str(batch_size),
        "-device", device,
        "-out_fpath", os.path.abspath(output_csv),
    ]
    proc = subprocess.run(
        cmd, cwd=str(CATAPRO_INFERENCE_DIR), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"CataPro predict.py failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        )
    if not os.path.isfile(output_csv):
        raise RuntimeError(
            f"CataPro produced no output file {output_csv}\n{proc.stdout}\n{proc.stderr}"
        )


def _reshape_native_output(native_csv: str, substrate: Optional[str]) -> pd.DataFrame:
    """Reshape CataPro's native (log10) CSV to the tps_eval schema keyed by ID."""
    native = pd.read_csv(native_csv)
    recovered_ids = [
        fasta_id[: -len(_ID_SUFFIX)] if str(fasta_id).endswith(_ID_SUFFIX) else str(fasta_id)
        for fasta_id in native[NATIVE_ID]
    ]
    frame = pd.DataFrame(
        {
            "ID": recovered_ids,
            "catapro_kcat": [log10_to_absolute(v) for v in native[NATIVE_KCAT]],
            "catapro_km": [log10_to_absolute(v) for v in native[NATIVE_KM]],
            "catapro_kcat_km": [log10_to_absolute(v) for v in native[NATIVE_KCAT_KM]],
        }
    )
    frame["catapro_substrate"] = substrate
    return frame[COLUMNS]


def score_fasta(
    fasta_path: str,
    target_substrate: Optional[str],
    *,
    smiles: Optional[str] = None,
    model_dpath: Optional[str] = None,
    device: str = "cuda",
    batch_size: int = 64,
    out_suffix: str = DEFAULT_OUT_SUFFIX,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """Predict CataPro kinetics for every sequence in a FASTA against one substrate.

    Writes a CSV keyed by ID (columns in COLUMNS) as a sibling of the FASTA
    (<fasta_stem>_<out_suffix>.csv) unless ``save_path`` is given. The output is
    reindexed to the full FASTA ID set (dropped/failed sequences -> NaN row). When
    no SMILES can be resolved for the substrate, every row is NaN (CataPro is not
    invoked).
    """
    ids, sequences = _read_fasta_ids_and_sequences(fasta_path)
    if save_path is None:
        save_path = _default_save_path(fasta_path, out_suffix)
    if model_dpath is None:
        model_dpath = str(DEFAULT_MODEL_DPATH)

    resolved_smiles = resolve_smiles(target_substrate, smiles)
    substrate_label = (target_substrate.upper() if target_substrate else None)

    if not resolved_smiles:
        print(
            f"[catapro] no SMILES for substrate {target_substrate!r}; "
            f"emitting NaN rows for {len(ids)} sequence(s)."
        )
        result = _nan_frame(ids, substrate_label)
        result.to_csv(save_path, index=False)
        print(f"[catapro] wrote {len(result)} rows to {save_path}")
        return result

    with tempfile.TemporaryDirectory(prefix="catapro_") as tmp:
        input_csv = os.path.join(tmp, "catapro_input.csv")
        native_csv = os.path.join(tmp, "catapro_native.csv")
        pd.DataFrame(
            {
                "Enzyme_id": ids,
                "type": CATAPRO_TYPE,
                "sequence": sequences,
                "smiles": resolved_smiles,
            }
        ).to_csv(input_csv)  # index written -> predict.py reads with index_col=0
        run_catapro_predict(
            input_csv,
            native_csv,
            model_dpath=model_dpath,
            device=device,
            batch_size=batch_size,
        )
        reshaped = _reshape_native_output(native_csv, substrate_label)

    # Reindex to the full FASTA ID set so any dropped sequence gets a NaN row.
    result = reshaped.set_index("ID").reindex(ids)
    result["catapro_substrate"] = substrate_label
    result = result.reset_index()[COLUMNS].sort_values("ID").reset_index(drop=True)
    result.to_csv(save_path, index=False)
    print(f"[catapro] wrote {len(result)} rows to {save_path} (substrate={substrate_label})")
    return result
