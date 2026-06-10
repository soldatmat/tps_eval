from __future__ import annotations

# ProteinMPNN sequence likelihood (NLL) — structure-branch metric.
#
# For each design structure, score its OWN sequence given its backbone with
# ProteinMPNN (`protein_mpnn_run.py --score_only 1`, no --path_to_fasta -> scores
# the native sequence read from the PDB). ProteinMPNN reports two negative
# log-likelihoods (averaged over residues):
#   * score        — over the DESIGNED residues only (here = all residues)
#   * global_score — over ALL residues in all chains
# For a single-chain monomer with everything designed these coincide; we report
# `proteinmpnn_nll` = global_score (the mean per-residue NLL). LOWER = the
# sequence is more compatible with / more likely given the fold.
#
# Implementation: shell out to the vendored ProteinMPNN (vendor/ProteinMPNN) once
# per structure with --score_only 1, which writes score_only/<stem>_pdb.npz with
# arrays `score` and `global_score` (shape (NUM_BATCHES,)). We average and read
# back. ID = structure filename stem (matches plddt / structural_identity keys).

import glob
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SRC_DIR.parent
PROTEINMPNN_DIR = REPO_ROOT / "vendor" / "ProteinMPNN"

COLUMNS = ["ID", "proteinmpnn_nll", "proteinmpnn_score_designed"]


def _collect_structures(structs_dir: str):
    """Map ID -> structure file. Mirrors plddt.py's af3-vs-flat detection.

    * "af3": AlphaFold3 af_output dir with per-job <job>/<job>_model.cif subfolders
      (ID = job name). ProteinMPNN reads CIF via parse_PDB (handles .cif).
    * "flat": directory of .pdb/.cif files (ID = filename stem); .pdb wins over .cif.
    """
    af3: Dict[str, str] = {}
    try:
        entries = sorted(os.listdir(structs_dir))
    except FileNotFoundError:
        entries = []
    for entry in entries:
        sub = os.path.join(structs_dir, entry)
        model = os.path.join(sub, entry + "_model.cif")
        if os.path.isdir(sub) and os.path.isfile(model):
            af3[entry] = model
    if af3:
        return OrderedDict(sorted(af3.items())), "af3"

    chosen: Dict[str, str] = {}
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, f"*{ext}"))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return OrderedDict(sorted(chosen.items())), "flat"


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_proteinmpnn_score.csv")


def score_pdb(
    pdb_path: str,
    out_folder: str,
    *,
    model_name: str = "v_48_020",
    seed: int = 0,
    num_passes: int = 1,
    backbone_noise: float = 0.0,
    python_exe: str = sys.executable,
) -> Dict[str, float]:
    """Run ProteinMPNN --score_only on a single structure's native sequence.
    Returns {'proteinmpnn_nll', 'proteinmpnn_score_designed'}."""
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    cmd = [
        python_exe,
        str(PROTEINMPNN_DIR / "protein_mpnn_run.py"),
        "--pdb_path", pdb_path,
        "--out_folder", out_folder,
        "--score_only", "1",
        "--num_seq_per_target", str(num_passes),
        "--seed", str(seed),
        "--batch_size", "1",
        "--backbone_noise", str(backbone_noise),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"ProteinMPNN failed on {pdb_path}:\n{proc.stdout}\n{proc.stderr}"
        )
    npz_path = os.path.join(out_folder, "score_only", stem + "_pdb.npz")
    if not os.path.isfile(npz_path):
        raise RuntimeError(
            f"ProteinMPNN produced no score file {npz_path}\n{proc.stdout}\n{proc.stderr}"
        )
    data = np.load(npz_path)
    return {
        "proteinmpnn_nll": float(np.mean(data["global_score"])),
        "proteinmpnn_score_designed": float(np.mean(data["score"])),
    }


def score_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    model_name: str = "v_48_020",
    seed: int = 0,
    backbone_noise: float = 0.0,
) -> pd.DataFrame:
    """Score every structure's native sequence with ProteinMPNN, writing a CSV
    keyed by ID (proteinmpnn_nll = mean per-residue NLL; lower = fold-compatible)."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected AF3 af_output or flat .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    with tempfile.TemporaryDirectory(prefix="proteinmpnn_score_") as tmp:
        for i, (stem, path) in enumerate(structures.items(), start=1):
            try:
                stats = score_pdb(
                    path, tmp, model_name=model_name, seed=seed, backbone_noise=backbone_noise
                )
            except Exception as exc:  # noqa: BLE001
                print(f"  [warn] failed on {stem}: {exc}")
                stats = {"proteinmpnn_nll": np.nan, "proteinmpnn_score_designed": np.nan}
            stats["ID"] = stem
            rows.append(stats)
            print(f"  [{i}/{n}] {stem}: nll={stats['proteinmpnn_nll']:.4f}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
