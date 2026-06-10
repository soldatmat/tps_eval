# -*- coding: utf-8 -*-
"""Aggrescan3D (A3D) structure-based aggregation-propensity metric for tps_eval.

A3D scores spatially-clustered, surface-exposed hydrophobic patches from a 3D
structure -- a structure-based expressibility/aggregation signal orthogonal to
the sequence-based SoluProt. Per residue, A3D emits an aggregation-propensity
score (positive = aggregation-prone, surface-exposed hydrophobic; negative =
protective/buried/polar). We run A3D in STATIC mode only (the dynamic CABS-flex
mode is far too slow for triage and is never triggered here -- we never pass
``--dynamic``) and reduce the per-residue A3D.csv to per-ID scalars.

IMPORTANT: this module is imported and run inside the dedicated ``aggrescan3d``
conda env, which is PYTHON 2.7 (Aggrescan3D's vendored source is Py2-only).
Keep this file Python-2.7 compatible: no f-strings, no py3-only typing, no
``from __future__ import annotations``.

Output: one row per input structure, keyed by ``ID`` (file stem), written to a
CSV that is a sibling of the structs dir, named ``<structs_dir>_aggregation.csv``.
Columns: ID, a3d_avg_score, a3d_total_score, a3d_max_score, a3d_min_score,
a3d_total_pos_score, n_residues. On A3D failure for a structure the row is
emitted with NaNs and processing continues.
"""

import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

# Output column order (ID first, then the filtration metrics).
COLUMNS = [
    "ID",
    "a3d_avg_score",
    "a3d_total_score",
    "a3d_max_score",
    "a3d_min_score",
    "a3d_total_pos_score",
    "n_residues",
]


def _is_cif(path):
    return path.lower().endswith((".cif", ".mmcif"))


def _collect_structures(structs_dir):
    """Map ID -> structure file. Mirrors plddt.py's flat/af3 auto-detection.

    * "af3"  -- an AlphaFold3 ``af_output`` directory: one subfolder per job,
      each with ``<job>/<job>_model.cif``. ID = job subfolder name.
    * "flat" -- a directory of ``.pdb``/``.cif`` files; ID = filename stem. If
      both a ``.pdb`` and ``.cif`` exist for an ID the ``.pdb`` wins (A3D wants
      PDB, so we avoid an unnecessary cif->pdb conversion).

    Returns (OrderedDict-like dict[id -> path] as sorted list of pairs, mode).
    """
    af3 = {}
    try:
        entries = sorted(os.listdir(structs_dir))
    except OSError:
        entries = []
    for entry in entries:
        sub = os.path.join(structs_dir, entry)
        model = os.path.join(sub, entry + "_model.cif")
        if os.path.isdir(sub) and os.path.isfile(model):
            af3[entry] = model
    if af3:
        return sorted(af3.items()), "af3"

    # Flat dir. Order matters: later wins, so .pdb is applied last.
    chosen = {}
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, "*" + ext))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return sorted(chosen.items()), "flat"


def _prepare_pdb(src_path, out_pdb):
    """Normalize any input structure (.pdb or .cif) into a clean PDB for A3D,
    using Biopython. Two things matter:

    1. A3D wants a PDB file; AF3 (and some flat) structures are mmCIF.
    2. Terminal ``OXT`` atoms break A3D's static scoring: the bundled freesasa
       (run with ``--radii naccess``) has no naccess radius for OXT, so the
       C-terminal residue's total SASA comes out ``-nan`` -> its relative SASA
       is printed as ``N/A`` -> A3D's ``float(...)`` on that column raises
       ``could not convert string to float: N/A`` and the whole structure
       fails. AlphaFold writes OXT on the C-terminus, so this hit real designs.
       OXT is just a terminal carboxyl oxygen and is irrelevant to the
       aggregation chemistry, so we drop it. (Observed on Aurum 2026-06: e.g.
       seq10114/seq1021 crashed, OXT strip makes them score fine.)

    Only the first model is kept (AF writes one; PDBIO writes a single MODEL).
    """
    from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    parser = MMCIFParser(QUIET=True) if _is_cif(src_path) else PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", src_path)
    # Keep only the first model so PDBIO writes a single MODEL.
    models = list(structure)
    for m in models[1:]:
        structure.detach_child(m.id)

    class _DropOXT(Select):
        def accept_atom(self, atom):
            return atom.get_name() != "OXT"

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb, select=_DropOXT())
    return out_pdb


def _parse_a3d_csv(csv_path):
    """Read A3D.csv (cols: protein,chain,residue,residue_name,score) -> np.array
    of per-residue scores (float). Empty array if no data rows."""
    scores = []
    with open(csv_path, "r") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            try:
                scores.append(float(row[4]))
            except (ValueError, IndexError):
                continue
    return np.asarray(scores, dtype=float)


def _summarize(scores):
    """Reduce a per-residue A3D score array to the per-ID scalar metrics."""
    arr = np.asarray(scores, dtype=float)
    if arr.size == 0:
        return {
            "a3d_avg_score": float("nan"),
            "a3d_total_score": float("nan"),
            "a3d_max_score": float("nan"),
            "a3d_min_score": float("nan"),
            "a3d_total_pos_score": float("nan"),
            "n_residues": 0,
        }
    return {
        "a3d_avg_score": float(arr.mean()),
        "a3d_total_score": float(arr.sum()),
        "a3d_max_score": float(arr.max()),
        "a3d_min_score": float(arr.min()),
        # Sum of positive contributions = total aggregation propensity (the
        # negative/protective residues don't reduce aggregation risk in vivo).
        "a3d_total_pos_score": float(arr[arr > 0].sum()),
        "n_residues": int(arr.size),
    }


def _run_a3d_static(pdb_path, work_dir):
    """Run Aggrescan3D in STATIC mode on a single PDB into work_dir; return the
    path to the produced A3D.csv. Raises CalledProcessError/RuntimeError on
    failure. STATIC mode = default (no --dynamic flag is ever passed)."""
    # -i input pdb, -w work_dir, -v 2 (quiet-ish; <4 so A3D cleans its own tmp).
    # No -d / --dynamic  => static SASA-based scoring with the bundled freesasa.
    cmd = [
        "aggrescan",
        "-i", pdb_path,
        "-w", work_dir,
        "-v", "2",
        "--overwrite",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _ = proc.communicate()
    if isinstance(out, bytes):
        out = out.decode("utf-8", "replace")
    csv_path = os.path.join(work_dir, "A3D.csv")
    if proc.returncode != 0 or not os.path.isfile(csv_path):
        raise RuntimeError(
            "aggrescan exited %s / no A3D.csv. tail:\n%s"
            % (proc.returncode, "\n".join(out.splitlines()[-15:]))
        )
    return csv_path


def _default_save_path(structs_dir):
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_aggregation.csv")


def extract_aggregation_dir(structs_dir, save_path=None, save_residue_scores=False,
                            residue_scores_dir=None):
    """Run A3D (static mode) on every structure in `structs_dir` and write a CSV
    keyed by ID with per-structure aggregation-propensity scalars.

    save_residue_scores: if True, also dump each structure's per-residue A3D
    score array to `residue_scores_dir` (one ``<ID>.csv`` per structure) for
    hotspot visualization. Default off.
    """
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            "No structures found in %s (expected an AlphaFold3 af_output dir with "
            "<job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif files)."
            % structs_dir
        )
    print("Detected %s layout: %d structure(s) in %s" % (mode, len(structures), structs_dir))

    if save_residue_scores:
        if residue_scores_dir is None:
            d = structs_dir.rstrip(os.sep)
            residue_scores_dir = os.path.join(
                os.path.dirname(d), os.path.basename(d) + "_aggregation_residue_scores"
            )
        if not os.path.isdir(residue_scores_dir):
            os.makedirs(residue_scores_dir)
        print("Per-residue scores -> %s" % residue_scores_dir)

    rows = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures, start=1):
        tmp_root = tempfile.mkdtemp(prefix="a3d_")
        try:
            # Always normalize through Biopython (handles cif->pdb AND strips the
            # OXT atoms that otherwise make freesasa emit N/A and crash A3D).
            pdb_path = os.path.join(tmp_root, stem + ".pdb")
            _prepare_pdb(path, pdb_path)
            work_dir = os.path.join(tmp_root, "run")
            csv_path = _run_a3d_static(pdb_path, work_dir)
            scores = _parse_a3d_csv(csv_path)
            stats = _summarize(scores)
            if save_residue_scores and scores.size:
                np.savetxt(
                    os.path.join(residue_scores_dir, stem + ".csv"),
                    scores, fmt="%.4f", delimiter=",",
                )
        except Exception as exc:  # A3D failure -> NaN row, keep going
            print("  [warn] A3D failed for %s: %s" % (os.path.basename(path), exc))
            stats = _summarize(np.asarray([], dtype=float))
            n_failed += 1
        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)
        stats["ID"] = stem
        rows.append(stats)
        print("  processed %d/%d (%s)" % (i, n, stem))
        sys.stdout.flush()

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print("Wrote %d rows to %s%s" % (
        len(df), save_path, (" (%d failed)" % n_failed) if n_failed else ""))
    return df
