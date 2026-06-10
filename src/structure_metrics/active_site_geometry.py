from __future__ import annotations

"""Active-site-specific structural metrics for class-I TPS designs.

Where ``motif_structural_distance`` measures the CA-CA *span* between the two
metal-binding motifs, this module looks at the side-chain-level geometry of the
catalytic carboxylate cage — i.e. whether the design's metal-coordinating oxygens
actually converge on a single locus that could hold the trinuclear Mg2+/Mn2+
cluster. It is fold-agnostic (AlphaFold or ESMFold .pdb/.cif) and apo-robust (no
metals or substrate need be modelled).

For each structure we

1. parse it with Biopython (reusing the af_output/flat auto-detection and the
   ID-stem convention from the neighbouring structure tools),
2. derive the 1-letter sequence + per-residue atom access,
3. run the SHARED motif localizer (``sequence_metrics.motif_localization``) to
   find the DDXXD-family and NSE/DTE motifs and their metal-coordinating residue
   offsets, and
4. gather the SIDE-CHAIN carboxylate/hydroxyl oxygens of those coordinating
   residues and compute the metrics below.

Metrics (columns, keyed by ``ID``):

* ``carboxylate_convergence_radius`` — RMS distance of the coordinating oxygens
  from their centroid. A competent trinuclear site clusters them within ~6-9 A;
  a splayed cage gives a large radius.
* ``n_coordinating_oxygens`` — count of carboxylate/hydroxyl O atoms found across
  both motifs (Asp OD1/OD2, Glu OE1/OE2, Asn OD1, Ser OG, Thr OG1).
* ``metal_point_void`` — clearance (A) at the oxygen centroid, i.e. the distance
  to the nearest PROTEIN atom *other than* the coordinating oxygens themselves. A
  real site leaves room (>~1.8 A) for the metals; a clash/filled centroid is bad.
* ``n_residues`` — modelled residue count, for context.

Geometric metrics are NaN when either motif is absent or no coordinating oxygen
is found (mirrors ``motif_structural_distance``'s NaN-on-absent-motif contract).

Optionally (``--templates ID[,ID...]``) we also report the **catalytic-
constellation RMSD**: a small 3D template of the catalytic constellation (Cα+Cβ
of the DDXXD + NSE/DTE coordinating residues, in fixed motif order) is built from
1-3 reference TPS structures in the same directory (e.g. ``1ps1``, ``5eat``); each
design's matched constellation is superposed onto every template and we report the
best (lowest) RMSD plus the winning template ID:

* ``catalytic_constellation_rmsd`` — best superposition RMSD (A) over the
  templates; NaN when motifs absent or a usable constellation can't be built.
* ``best_template`` — ID of the template that gave the best RMSD ('' when NaN).

NOTE: PyMOL's ``super``/``align`` would be the preferred superposer, but PyMOL is
not importable in the tps_eval env, so we use Biopython's ``Superimposer`` on the
matched Cα+Cβ atoms (a fixed-cardinality, ordered correspondence) instead.
"""

import argparse
import glob
import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Reuse the sequence-metrics shared motif localization (single source of truth).
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from sequence_metrics.motif_localization import (  # noqa: E402
    DDXXD_COORDINATING_OFFSETS,
    NSE_DTE_COORDINATING_OFFSETS,
    coordinating_indices,
    locate_ddxxd,
    locate_nse_dte,
)

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

COLUMNS = [
    "ID",
    "carboxylate_convergence_radius",
    "n_coordinating_oxygens",
    "metal_point_void",
    "catalytic_constellation_rmsd",
    "best_template",
    "n_residues",
]

# Backbone atoms used to define the catalytic-constellation correspondence. Cα+Cβ
# give a fixed-cardinality, ordered atom set per coordinating residue (glycine has
# no CB, but the coordinating residues are never Gly — they're D/E/N/S/T).
CONSTELLATION_ATOMS: Tuple[str, ...] = ("CA", "CB")

# Side-chain carboxylate / hydroxyl oxygen atoms that coordinate the metal cluster,
# per residue type (Biopython atom names). Only these atoms are gathered for the
# convergence-radius / centroid computation.
COORDINATING_OXYGEN_ATOMS: Dict[str, Tuple[str, ...]] = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "ASN": ("OD1",),
    "SER": ("OG",),
    "THR": ("OG1",),
}


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def _three_to_one(resname: str) -> str:
    """3-letter residue name -> 1-letter code, 'X' for anything non-standard."""
    try:
        return index_to_one(three_to_index(resname))
    except KeyError:
        return "X"


def structure_sequence_residues_atoms(
    structure_path: str,
) -> Tuple[str, list, np.ndarray]:
    """Parse one structure and return ``(sequence, residues, all_atom_coords)``.

    * ``sequence`` — 1-letter sequence of the protein residues of the first model,
      in chain order (HETATM ions/ligands/water skipped), matching the other
      structure tools' sequence derivation.
    * ``residues`` — the index-aligned list of Biopython ``Residue`` objects (so
      we can reach into their side-chain atoms by name).
    * ``all_atom_coords`` — (N, 3) array of every PROTEIN atom coordinate in the
      first model (used for the metal-point-void clearance check).
    """
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))  # first model only (predicted structs write one)

    seq_chars: List[str] = []
    residues: list = []
    all_atom_coords: List[np.ndarray] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATM (ions/ligands/water)
                continue
            seq_chars.append(_three_to_one(residue.get_resname()))
            residues.append(residue)
            for atom in residue:
                all_atom_coords.append(np.asarray(atom.get_coord(), dtype=float))
    coords = np.vstack(all_atom_coords) if all_atom_coords else np.empty((0, 3))
    return "".join(seq_chars), residues, coords


def _coordinating_oxygens(indices: List[int], residues: list) -> np.ndarray:
    """Gather side-chain carboxylate/hydroxyl O coords for the given residue
    indices. Residues whose expected atoms are absent (or that aren't an
    oxygen-bearing coordinating type) are skipped. Returns an (n, 3) array."""
    pts: List[np.ndarray] = []
    for i in indices:
        if not (0 <= i < len(residues)):
            continue
        residue = residues[i]
        for atom_name in COORDINATING_OXYGEN_ATOMS.get(residue.get_resname(), ()):
            if atom_name in residue:
                pts.append(np.asarray(residue[atom_name].get_coord(), dtype=float))
    return np.vstack(pts) if pts else np.empty((0, 3))


def _coordinating_indices_both(sequence: str) -> Optional[List[int]]:
    """0-based sequence indices of the metal-coordinating residues of BOTH motifs
    (DDXXD then NSE/DTE), or None if either motif is absent."""
    ddxxd = locate_ddxxd(sequence)
    nse = locate_nse_dte(sequence)
    if ddxxd is None or nse is None:
        return None
    return coordinating_indices(ddxxd, DDXXD_COORDINATING_OFFSETS) + coordinating_indices(
        nse, NSE_DTE_COORDINATING_OFFSETS
    )


def constellation_atoms(
    sequence: str, residues: list
) -> Optional["OrderedDict[Tuple[int, str], np.ndarray]"]:
    """Ordered map (residue-index, atom-name) -> coordinate for the Cα+Cβ atoms of
    the coordinating residues of both motifs, in fixed motif/atom order. The keys
    encode the *motif slot* (position within the concatenated coordinating list, not
    the absolute residue index) so a design and a template are matched slot-by-slot.
    Returns None when either motif is absent."""
    idx = _coordinating_indices_both(sequence)
    if idx is None:
        return None
    atoms: "OrderedDict[Tuple[int, str], np.ndarray]" = OrderedDict()
    for slot, i in enumerate(idx):
        if not (0 <= i < len(residues)):
            continue
        residue = residues[i]
        for atom_name in CONSTELLATION_ATOMS:
            if atom_name in residue:
                atoms[(slot, atom_name)] = np.asarray(residue[atom_name].get_coord(), dtype=float)
    return atoms


def _superpose_rmsd(
    template: "OrderedDict[Tuple[int, str], np.ndarray]",
    query: "OrderedDict[Tuple[int, str], np.ndarray]",
) -> Optional[float]:
    """RMSD of superposing ``query`` onto ``template`` over their shared
    (slot, atom) keys via Biopython's Superimposer. None if <3 atoms in common."""
    shared = [k for k in template if k in query]
    if len(shared) < 3:
        return None
    from Bio.PDB.Atom import Atom

    def _atoms(src):
        return [
            Atom(k[1], src[k], 0.0, 1.0, " ", k[1], i, element="C")
            for i, k in enumerate(shared)
        ]

    sup = Superimposer()
    sup.set_atoms(_atoms(template), _atoms(query))
    return float(sup.rms)


def build_templates(
    structs_dir: str, template_ids: List[str]
) -> "OrderedDict[str, OrderedDict[Tuple[int, str], np.ndarray]]":
    """Build catalytic-constellation templates for the given reference IDs from the
    structures in ``structs_dir``. IDs whose structure can't be parsed or whose
    motifs aren't found are skipped (with a warning)."""
    structures, _ = _collect_structures(structs_dir)
    templates: "OrderedDict[str, OrderedDict[Tuple[int, str], np.ndarray]]" = OrderedDict()
    for tid in template_ids:
        path = structures.get(tid)
        if path is None:
            print(f"  [warn] template '{tid}' not found in {structs_dir}; skipping")
            continue
        try:
            seq, residues, _ = structure_sequence_residues_atoms(path)
            atoms = constellation_atoms(seq, residues)
        except Exception as exc:
            print(f"  [warn] template '{tid}' failed to parse ({exc}); skipping")
            continue
        if not atoms:
            print(f"  [warn] template '{tid}' has no locatable catalytic motifs; skipping")
            continue
        templates[tid] = atoms
    return templates


def active_site_geometry(
    structure_path: str,
    templates: Optional["OrderedDict[str, OrderedDict[Tuple[int, str], np.ndarray]]"] = None,
) -> Dict[str, float]:
    """Carboxylate-cage geometry metrics for one structure; NaN geometry when a
    motif is absent or no coordinating oxygen is found. If ``templates`` is given,
    also reports the best catalytic-constellation RMSD and the winning template."""
    sequence, residues, all_atom_coords = structure_sequence_residues_atoms(structure_path)
    result: Dict[str, float] = {
        "carboxylate_convergence_radius": np.nan,
        "n_coordinating_oxygens": 0,
        "metal_point_void": np.nan,
        "catalytic_constellation_rmsd": np.nan,
        "best_template": "",
        "n_residues": len(sequence),
    }

    ddxxd = locate_ddxxd(sequence)
    nse = locate_nse_dte(sequence)
    if ddxxd is None or nse is None:
        return result

    idx = coordinating_indices(ddxxd, DDXXD_COORDINATING_OFFSETS) + coordinating_indices(
        nse, NSE_DTE_COORDINATING_OFFSETS
    )
    oxygens = _coordinating_oxygens(idx, residues)
    result["n_coordinating_oxygens"] = int(len(oxygens))
    if len(oxygens) == 0:
        return result

    centroid = oxygens.mean(axis=0)
    # RMS distance of the coordinating oxygens from their centroid.
    radius = float(np.sqrt(((oxygens - centroid) ** 2).sum(axis=1).mean()))
    result["carboxylate_convergence_radius"] = radius

    # Metal-point void: clearance at the centroid to the nearest protein atom that
    # is NOT one of the coordinating oxygens. Exclude the coordinating oxygens by
    # coordinate-matching (they are a subset of all_atom_coords).
    if len(all_atom_coords):
        d = np.sqrt(((all_atom_coords - centroid) ** 2).sum(axis=1))
        # Mask out the coordinating oxygens themselves (exact coord match).
        keep = np.ones(len(all_atom_coords), dtype=bool)
        for ox in oxygens:
            keep &= ~np.all(np.isclose(all_atom_coords, ox), axis=1)
        if keep.any():
            result["metal_point_void"] = float(d[keep].min())

    # Catalytic-constellation RMSD against the reference templates (Biopython
    # Superimposer fallback; PyMOL super/align unavailable in this env).
    if templates:
        query = constellation_atoms(sequence, residues)
        if query:
            best_rmsd = np.inf
            best_tid = ""
            for tid, tmpl in templates.items():
                rmsd = _superpose_rmsd(tmpl, query)
                if rmsd is not None and rmsd < best_rmsd:
                    best_rmsd, best_tid = rmsd, tid
            if np.isfinite(best_rmsd):
                result["catalytic_constellation_rmsd"] = float(best_rmsd)
                result["best_template"] = best_tid
    return result


def _collect_structures(structs_dir: str) -> Tuple["OrderedDict[str, str]", str]:
    """Map ID -> structure file, auto-detecting layout. Mirrors
    motif_structural_distance/plddt: an AF3 ``af_output`` dir (per-job
    ``<job>/<job>_model.cif``; ID = job name) takes precedence; otherwise a flat
    dir of .pdb/.cif (ID = filename stem; .pdb wins on tie)."""
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
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_active_site_geometry.csv")


def active_site_geometry_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    template_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Carboxylate-cage geometry for every structure in ``structs_dir``; CSV keyed by ID.

    If ``template_ids`` is given, also computes the catalytic-constellation RMSD of
    every structure against templates built from those reference IDs."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected an AlphaFold3 af_output "
            "dir with <job>/<job>_model.cif subfolders, or a flat dir of .pdb/.cif)."
        )
    print(f"Detected {mode} layout: {len(structures)} structure(s) in {structs_dir}")

    templates = None
    if template_ids:
        templates = build_templates(structs_dir, template_ids)
        if templates:
            print(f"Built {len(templates)} catalytic-constellation template(s): "
                  f"{', '.join(templates)}")
        else:
            print("[warn] no usable templates built; catalytic_constellation_rmsd will be NaN")

    rows: List[Dict[str, float]] = []
    n = len(structures)
    n_failed = 0
    for i, (stem, path) in enumerate(structures.items(), start=1):
        try:
            stats = active_site_geometry(path, templates=templates)
        except Exception as exc:  # malformed/unparsable -> NaN row, keep going
            print(f"  [warn] failed to parse {os.path.basename(path)}: {exc}")
            stats = {
                "carboxylate_convergence_radius": np.nan,
                "n_coordinating_oxygens": 0,
                "metal_point_void": np.nan,
                "catalytic_constellation_rmsd": np.nan,
                "best_template": "",
                "n_residues": 0,
            }
            n_failed += 1
        stats["ID"] = str(stem).strip()
        rows.append(stats)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}" + (f" ({n_failed} unparsable)" if n_failed else ""))
    return df
