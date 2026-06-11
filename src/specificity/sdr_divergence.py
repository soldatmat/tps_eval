from __future__ import annotations

"""Specificity-determining-residue (SDR) divergence for class-I TPS designs.

The failure mode this catches: a generated design is GLOBALLY close to a known
product-characterized TPS (high TM-score / sequence identity) yet DIVERGES at the
handful of active-site residues that determine *which* terpene the enzyme makes.
The textbook case is TEAS vs HPS — ~75-80% identical overall, indistinguishable by
global similarity, yet making different sesquiterpenes because of ~9 active-site
residues (Greenhagen/O'Maille, PNAS 2006). A k-NN / global-similarity label transfer
would happily copy the neighbour's product label onto such a design; this metric
flags it instead.

Per design we:

1. Take its NEAREST known-TPS neighbour + global similarity from the committed
   ``--top_k`` outputs (``structural_identity`` TM-score and/or
   ``max_sequence_identity`` percent), rank-1. The structural ``neighbour_id`` may
   carry a foldseek ``_<chain>`` suffix — it is stripped to map back to the
   reference structure/sequence (reusing ``knn``'s ``_strip_chain_suffix``).
2. Locate the SDR positions on that neighbour (see SDR PANEL below) and
   STRUCTURALLY ALIGN the design onto the neighbour (Biopython ``Superimposer`` on
   anchor Cα atoms, the same env-safe approach ``active_site_geometry`` uses — PyMOL
   is unavailable). Each neighbour SDR residue is mapped to the design residue whose
   Cα is nearest after superposition (within a tolerance).
3. Compute ``sdr_identity`` = fraction of mapped SDR positions where the design
   residue equals the neighbour residue; record ``n_sdr_positions``,
   ``n_sdr_mismatches``, the compact ``divergent_positions`` string, and the flag
   ``specificity_divergence`` = (global similarity >= tau_high) AND
   (sdr_identity <= tau_low) — "looks like this known TPS overall, but differs where
   it matters".

SDR PANEL — modular (an input), like the k-NN's label file:

* Default (structure-derived, no panel file): the active-site-lining residues, i.e.
  every residue with a heavy atom within ``--panel_cutoff`` A (default 10 A) of the
  carboxylate-cage metal point computed from the neighbour's own structure (the
  centroid of the DDXXD coordinating oxygens, plus the NSE/DTE ones when that motif is
  also matched — see ``metal_point`` for why NSE/DTE is optional here). This is close
  to the active-site residue set the structure-based
  classifier uses; it needs no literature position-mapping and adapts per neighbour.
* Explicit override (``--sdr_panel <file>``): a committable CSV of explicit SDR
  positions anchored to a named reference (see ``load_panel`` for the format). The
  tool is agnostic to its contents.

Output ``<structs_dir>_sdr_divergence.csv`` keyed by ``ID``:
    nearest_neighbour_id, nearest_neighbour_similarity, similarity_space,
    n_sdr_positions, sdr_identity, n_sdr_mismatches,
    specificity_divergence (bool), divergent_positions

NaN / flag=False when no neighbour clears tau_high, or when the SDR positions can't
be located/mapped (missing motifs for the structure-derived panel, unparsable
structure, or fewer than 3 superposition anchors).
"""

import glob
import math
import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import index_to_one, three_to_index
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Reuse the sibling structure-metrics + sequence-metrics modules (single source of
# truth for the metal-point geometry, motif localization, and chain-suffix strip).
SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from structure_metrics.active_site_geometry import (  # noqa: E402
    metal_point as _cage_metal_point,
)
try:
    # Reuse the k-NN tool's chain-suffix strip when it's present (single source of
    # truth). It is built in parallel and may not be checked out yet, so fall back to
    # an inline copy with the identical contract.
    from knn.knn_label_transfer import _strip_chain_suffix  # noqa: E402
except Exception:  # pragma: no cover - exercised only when knn is absent
    def _strip_chain_suffix(neighbour_id, valid_ids=None):  # type: ignore
        """Strip a foldseek ``_<chain>`` suffix from a structural neighbour id.

        Mirror of ``knn.knn_label_transfer._strip_chain_suffix``: with ``valid_ids``,
        strip a trailing ``_<token>`` only when the full id is unknown but the stem is
        known; without it, strip a short (<=2 char) alphanumeric trailing token."""
        nid = str(neighbour_id)
        if valid_ids is not None:
            if nid in valid_ids:
                return nid
            if "_" in nid:
                stem = nid.rsplit("_", 1)[0]
                if stem in valid_ids:
                    return stem
            return nid
        if "_" in nid:
            stem, tok = nid.rsplit("_", 1)
            if len(tok) <= 2 and tok.isalnum():
                return stem
        return nid

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)

# --------------------------------------------------------------------------- #
# tau defaults (documented; grounded in the TEAS/HPS literature)
# --------------------------------------------------------------------------- #
# TEAS and HPS are ~75-80% identical yet differ at ~9 active-site residues and make
# different products (Greenhagen et al., PNAS 2006). So the regime we want to flag is
# "globally very similar but locally divergent": a permissive global floor (well below
# the 75-80% they share, so any near-neighbour of that calibre clears it) combined
# with a strict local-identity ceiling (an active-site panel sharing <= ~70% of its
# residues with the neighbour is a meaningful specificity divergence).
#   tau_high : global-similarity floor for "looks like this known TPS overall".
#              On the SIMILARITY scale in [0,1] (TM-score, or identity%/100).
#              0.6 ~= TM 0.6 / 60% id  (comfortably above the fold-similarity floor
#              of TM 0.5, and below the TEAS/HPS shared ~0.78).
#   tau_low  : SDR-identity ceiling for "differs where it matters". 0.7 -> a design
#              matching <= 70% of the panel diverges at >= 30% of specificity residues.
DEFAULT_TAU_HIGH = 0.6
DEFAULT_TAU_LOW = 0.7

# Default structure-derived panel cutoff (A) around the metal point. 8-10 A captures
# the active-site-lining first/second shell; 10 A is the inclusive default.
DEFAULT_PANEL_CUTOFF = 10.0

# Tolerance (A) for accepting a design residue as the structural counterpart of a
# neighbour SDR residue after Cα superposition. Generous enough to absorb backbone
# drift between a design and its homolog, tight enough to reject a non-correspondence.
DEFAULT_MAP_TOLERANCE = 4.0

COLUMNS = [
    "ID",
    "nearest_neighbour_id",
    "nearest_neighbour_similarity",
    "similarity_space",
    "n_sdr_positions",
    "sdr_identity",
    "n_sdr_mismatches",
    "specificity_divergence",
    "divergent_positions",
]


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def _three_to_one(resname: str) -> str:
    try:
        return index_to_one(three_to_index(resname))
    except KeyError:
        return "X"


# --------------------------------------------------------------------------- #
# Structure collection (mirrors plddt/active_site_geometry)
# --------------------------------------------------------------------------- #
def collect_structures(structs_dir: str) -> Tuple["OrderedDict[str, str]", str]:
    """Map ID -> structure file, auto-detecting AF3 ``af_output`` vs a flat dir.

    Identical contract to ``active_site_geometry._collect_structures`` /
    ``plddt``: an AF3 ``af_output`` dir (per-job ``<job>/<job>_model.cif``; ID = job
    name) takes precedence; otherwise a flat dir of .pdb/.cif (ID = filename stem;
    .pdb wins on tie)."""
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


# --------------------------------------------------------------------------- #
# Per-residue structural representation
# --------------------------------------------------------------------------- #
class ResidueInfo:
    """Compact per-residue view of a parsed structure for SDR work.

    Index-aligned lists (by 0-based protein-residue order, the same order
    ``structure_sequence_residues_atoms`` returns):
      seq          1-letter sequence (HETATM skipped)
      residues     Biopython Residue objects
      ca           (N,3) Cα coords (NaN row if a residue has no CA)
      heavy        list of (n_i,3) arrays of that residue's heavy-atom coords
      resnums      author residue numbers (int; for panel anchoring / reporting)
    """

    def __init__(self, structure_path: str):
        parser = _parser_for(structure_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure("s", structure_path)
        model = next(iter(structure))
        seq_chars: List[str] = []
        residues: list = []
        ca: List[np.ndarray] = []
        heavy: List[np.ndarray] = []
        resnums: List[int] = []
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                seq_chars.append(_three_to_one(residue.get_resname()))
                residues.append(residue)
                resnums.append(int(residue.id[1]))
                if "CA" in residue:
                    ca.append(np.asarray(residue["CA"].get_coord(), dtype=float))
                else:
                    ca.append(np.full(3, np.nan))
                pts = [
                    np.asarray(a.get_coord(), dtype=float)
                    for a in residue
                    if a.element != "H"
                ]
                heavy.append(np.vstack(pts) if pts else np.empty((0, 3)))
        self.seq = "".join(seq_chars)
        self.residues = residues
        self.ca = np.vstack(ca) if ca else np.empty((0, 3))
        self.heavy = heavy
        self.resnums = resnums


# --------------------------------------------------------------------------- #
# Metal point + structure-derived panel
# --------------------------------------------------------------------------- #
def metal_point(info: ResidueInfo) -> Optional[np.ndarray]:
    """Carboxylate-cage metal point for the SDR panel anchor — a thin ``ResidueInfo``
    adapter over the CANONICAL ``active_site_geometry.metal_point`` (the single source
    of truth). Centroid of the DDXXD (+ NSE/DTE when that motif is also matched)
    coordinating side-chain oxygens.

    DDXXD is required and NSE/DTE is an optional refinement: the DDXXD coordinating
    oxygens alone localize the catalytic centre well, and the shared NSE/DTE regex does
    NOT match real TEAS (5EAT), so requiring it would make the structure-derived panel
    fail on the very reference TPS we target. None only if DDXXD is absent or no DDXXD
    coordinating oxygen is found."""
    return _cage_metal_point(info.seq, info.residues)


def structure_derived_panel(info: ResidueInfo, cutoff: float) -> Optional[List[int]]:
    """0-based residue indices of the active-site-lining residues: those with any
    heavy atom within ``cutoff`` A of the neighbour's metal point. None if the metal
    point can't be computed (motifs absent)."""
    mp = metal_point(info)
    if mp is None:
        return None
    idx: List[int] = []
    c2 = cutoff * cutoff
    for i, hv in enumerate(info.heavy):
        if hv.shape[0] == 0:
            continue
        d2 = ((hv - mp) ** 2).sum(axis=1)
        if (d2 <= c2).any():
            idx.append(i)
    return idx


# --------------------------------------------------------------------------- #
# Explicit SDR panel file
# --------------------------------------------------------------------------- #
def load_panel(panel_file: str) -> "Panel":
    """Load an explicit SDR panel CSV anchored to a named reference.

    Format (header required), one row per SDR position::

        reference_id,resnum,expected_residue
        TEAS_5eat,274,T
        TEAS_5eat,294,T
        ...

    Columns:
      reference_id     The reference structure/sequence the resnums are anchored to.
                       A panel may cover several references (e.g. one per known class);
                       only rows whose reference_id matches a design's nearest
                       neighbour id (after chain-suffix strip) are used for it.
      resnum           Author residue number in that reference structure.
      expected_residue (optional) 1-letter code expected at that position; informational
                       only — the metric compares design-vs-neighbour, not vs expected.

    A panel with NO row matching a given neighbour id contributes no positions for that
    design (it falls back to NaN, NOT to the structure-derived default — mixing the two
    silently would be misleading; choose one mode per run).

    Lines beginning with ``#`` are treated as comments (so a committed panel can carry
    inline provenance/caveats); any extra columns beyond the recognized ones are ignored.
    """
    df = pd.read_csv(panel_file, comment="#", skip_blank_lines=True)
    cols = {c.lower(): c for c in df.columns}
    if "reference_id" not in cols or "resnum" not in cols:
        raise ValueError(
            f"SDR panel {panel_file} must have at least 'reference_id' and 'resnum' columns; "
            f"got {list(df.columns)}"
        )
    ref_col, num_col = cols["reference_id"], cols["resnum"]
    exp_col = cols.get("expected_residue")
    by_ref: Dict[str, List[Tuple[int, Optional[str]]]] = {}
    for _, row in df.iterrows():
        rid = str(row[ref_col]).strip()
        try:
            resnum = int(row[num_col])
        except (ValueError, TypeError):
            continue
        exp = None
        if exp_col is not None and not pd.isna(row[exp_col]):
            exp = str(row[exp_col]).strip()[:1].upper() or None
        by_ref.setdefault(rid, []).append((resnum, exp))
    return Panel(by_ref)


class Panel:
    """Explicit SDR panel: reference_id -> list of (author resnum, expected residue)."""

    def __init__(self, by_ref: Dict[str, List[Tuple[int, Optional[str]]]]):
        self.by_ref = by_ref

    def references(self) -> List[str]:
        return sorted(self.by_ref)

    def indices_for(self, neighbour_id: str, info: ResidueInfo) -> Optional[List[int]]:
        """0-based residue indices in ``info`` for the panel rows anchored to
        ``neighbour_id`` (matched against panel reference ids, allowing the panel id
        to be a prefix or the neighbour id to be a prefix, e.g. structure stem vs
        full id). None if the neighbour matches no panel reference; [] if it matches
        but none of its resnums are present in the structure."""
        ref = self._match_reference(neighbour_id)
        if ref is None:
            return None
        num_to_idx = {rn: i for i, rn in enumerate(info.resnums)}
        idx: List[int] = []
        for resnum, _exp in self.by_ref[ref]:
            i = num_to_idx.get(resnum)
            if i is not None:
                idx.append(i)
        return sorted(set(idx))

    def _match_reference(self, neighbour_id: str) -> Optional[str]:
        nid = str(neighbour_id)
        if nid in self.by_ref:
            return nid
        # Tolerate stem/full-id mismatches in either direction.
        for ref in self.by_ref:
            if nid == ref or nid.startswith(ref) or ref.startswith(nid):
                return ref
        return None


# --------------------------------------------------------------------------- #
# Superposition + position mapping
# --------------------------------------------------------------------------- #
def _seq_anchor_pairs(
    design: ResidueInfo, neighbour: ResidueInfo
) -> List[Tuple[int, int]]:
    """Anchor (design_idx, neighbour_idx) pairs for superposition via a global
    sequence alignment of the two sequences (Cα atoms of aligned, identical-or-not
    positions). Uses Biopython's PairwiseAligner; falls back to a positional 1:1 map
    on equal-length sequences if the aligner is unavailable."""
    ds, ns = design.seq, neighbour.seq
    try:
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aln = aligner.align(ds, ns)[0]
        pairs: List[Tuple[int, int]] = []
        for (d0, d1), (n0, n1) in zip(aln.aligned[0], aln.aligned[1]):
            for off in range(d1 - d0):
                pairs.append((d0 + off, n0 + off))
        return pairs
    except Exception:
        n = min(len(ds), len(ns))
        return [(i, i) for i in range(n)]


def _superpose(
    design: ResidueInfo, neighbour: ResidueInfo, pairs: Sequence[Tuple[int, int]]
) -> Optional[Superimposer]:
    """Biopython Superimposer fit of the design onto the neighbour over the Cα atoms
    of the aligned pairs (PyMOL super/align is unavailable in this env). None if
    fewer than 3 usable (non-NaN) Cα pairs."""
    d_atoms: List[Atom] = []
    n_atoms: List[Atom] = []
    for k, (di, ni) in enumerate(pairs):
        if not (0 <= di < len(design.ca) and 0 <= ni < len(neighbour.ca)):
            continue
        dc, nc = design.ca[di], neighbour.ca[ni]
        if np.isnan(dc).any() or np.isnan(nc).any():
            continue
        n_atoms.append(Atom("CA", nc, 0.0, 1.0, " ", "CA", k, element="C"))
        d_atoms.append(Atom("CA", dc, 0.0, 1.0, " ", "CA", k, element="C"))
    if len(d_atoms) < 3:
        return None
    sup = Superimposer()
    sup.set_atoms(n_atoms, d_atoms)  # fixed = neighbour, moving = design
    return sup


def _map_neighbour_to_design(
    sup: Superimposer,
    design: ResidueInfo,
    neighbour: ResidueInfo,
    neighbour_idx: Sequence[int],
    tolerance: float,
) -> Dict[int, int]:
    """For each neighbour SDR residue index, find the design residue whose Cα is
    nearest (after superposing the design onto the neighbour frame) and within
    ``tolerance`` A. Returns {neighbour_idx: design_idx} for the mapped ones."""
    rot, tran = sup.rotran
    # Move design Cα into the neighbour frame.
    d_ca = design.ca
    valid = ~np.isnan(d_ca).any(axis=1)
    moved = np.full_like(d_ca, np.nan)
    moved[valid] = d_ca[valid] @ rot + tran
    mapping: Dict[int, int] = {}
    for ni in neighbour_idx:
        if not (0 <= ni < len(neighbour.ca)):
            continue
        nc = neighbour.ca[ni]
        if np.isnan(nc).any():
            continue
        d2 = ((moved - nc) ** 2).sum(axis=1)
        d2[~valid] = np.inf
        j = int(np.argmin(d2))
        if math.sqrt(d2[j]) <= tolerance:
            mapping[ni] = j
    return mapping


# --------------------------------------------------------------------------- #
# Per-design SDR divergence
# --------------------------------------------------------------------------- #
def _to_similarity(space: str, score: float) -> float:
    """Global-similarity score -> [0,1]. structural=TM-score (as-is); sequence=%/100."""
    if score is None or (isinstance(score, float) and math.isnan(score)):
        return float("nan")
    return float(score) / 100.0 if space == "sequence" else float(score)


def sdr_divergence_one(
    design_info: ResidueInfo,
    neighbour_id: str,
    neighbour_info: ResidueInfo,
    similarity: float,
    space: str,
    *,
    panel: Optional[Panel] = None,
    panel_cutoff: float = DEFAULT_PANEL_CUTOFF,
    map_tolerance: float = DEFAULT_MAP_TOLERANCE,
    tau_high: float = DEFAULT_TAU_HIGH,
    tau_low: float = DEFAULT_TAU_LOW,
) -> Dict[str, object]:
    """SDR divergence of one design against its nearest neighbour.

    ``similarity`` is already on the [0,1] scale. Returns the row dict (sans ID)."""
    row: Dict[str, object] = {
        "nearest_neighbour_id": neighbour_id,
        "nearest_neighbour_similarity": float(similarity),
        "similarity_space": space,
        "n_sdr_positions": 0,
        "sdr_identity": np.nan,
        "n_sdr_mismatches": np.nan,
        "specificity_divergence": False,
        "divergent_positions": "",
    }

    # SDR positions on the NEIGHBOUR (reference).
    if panel is not None:
        nbr_idx = panel.indices_for(neighbour_id, neighbour_info)
    else:
        nbr_idx = structure_derived_panel(neighbour_info, panel_cutoff)
    if not nbr_idx:
        return row  # no panel positions -> can't map -> NaN

    # Structurally align design onto neighbour and map the SDR residues.
    pairs = _seq_anchor_pairs(design_info, neighbour_info)
    sup = _superpose(design_info, neighbour_info, pairs)
    if sup is None:
        return row
    mapping = _map_neighbour_to_design(
        sup, design_info, neighbour_info, nbr_idx, map_tolerance
    )
    if not mapping:
        return row

    n_pos = len(mapping)
    mism: List[str] = []
    n_match = 0
    for ni, dj in sorted(mapping.items()):
        nbr_res = neighbour_info.seq[ni] if 0 <= ni < len(neighbour_info.seq) else "X"
        des_res = design_info.seq[dj] if 0 <= dj < len(design_info.seq) else "X"
        resnum = neighbour_info.resnums[ni] if 0 <= ni < len(neighbour_info.resnums) else ni + 1
        if des_res == nbr_res:
            n_match += 1
        else:
            mism.append(f"{nbr_res}{resnum}{des_res}")

    sdr_identity = n_match / n_pos if n_pos else float("nan")
    row["n_sdr_positions"] = int(n_pos)
    row["sdr_identity"] = float(sdr_identity)
    row["n_sdr_mismatches"] = int(n_pos - n_match)
    row["divergent_positions"] = ";".join(mism)
    row["specificity_divergence"] = bool(
        (not math.isnan(similarity))
        and similarity >= tau_high
        and (not math.isnan(sdr_identity))
        and sdr_identity <= tau_low
    )
    return row


# --------------------------------------------------------------------------- #
# Nearest-neighbour selection from top-k CSVs
# --------------------------------------------------------------------------- #
def _rank1_neighbours(
    structural_topk: Optional[str],
    sequence_topk: Optional[str],
    valid_ref_ids: set,
) -> Dict[str, Tuple[str, float, str]]:
    """{design_id: (neighbour_id, similarity, space)} from the rank-1 row of the
    top-k CSVs. Structural is preferred when present (and clears NaN); falls back to
    sequence. The structural neighbour_id is chain-suffix-stripped against the known
    reference ids. Designs absent from both CSVs are omitted (caller NaNs them)."""
    best: Dict[str, Tuple[str, float, str]] = {}

    def _ingest(path: Optional[str], space: str, prefer: bool) -> None:
        if not path:
            return
        df = pd.read_csv(path)
        df = df.sort_values(["query_id", "rank"], kind="stable")
        for qid, grp in df.groupby("query_id", sort=False):
            top = grp.iloc[0]
            nid = str(top["neighbour_id"])
            if space == "structural":
                nid = _strip_chain_suffix(nid, valid_ref_ids)
            sim = _to_similarity(space, top["score"])
            qid = str(qid)
            if math.isnan(sim):
                continue
            if prefer or qid not in best:
                best[qid] = (nid, sim, space)

    # Sequence first (fills), then structural overrides (preferred).
    _ingest(sequence_topk, "sequence", prefer=False)
    _ingest(structural_topk, "structural", prefer=True)
    return best


# --------------------------------------------------------------------------- #
# Directory driver
# --------------------------------------------------------------------------- #
def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_sdr_divergence.csv")


def sdr_divergence_dir(
    structs_dir: str,
    known_structs_dir: str,
    *,
    structural_topk: Optional[str] = None,
    sequence_topk: Optional[str] = None,
    panel_file: Optional[str] = None,
    panel_cutoff: float = DEFAULT_PANEL_CUTOFF,
    map_tolerance: float = DEFAULT_MAP_TOLERANCE,
    tau_high: float = DEFAULT_TAU_HIGH,
    tau_low: float = DEFAULT_TAU_LOW,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """SDR divergence for every design in ``structs_dir`` against its nearest
    known-TPS neighbour (rank-1 from the top-k CSVs), written to a CSV keyed by ID.

    Requires at least one of ``structural_topk`` / ``sequence_topk``. ``known_structs_dir``
    supplies the reference structures the neighbour ids resolve to."""
    if not structural_topk and not sequence_topk:
        raise ValueError(
            "Provide at least one of structural_topk / sequence_topk to identify the "
            "nearest neighbour per design."
        )

    designs, dmode = collect_structures(structs_dir)
    if not designs:
        raise ValueError(f"No design structures found in {structs_dir}.")
    knowns, kmode = collect_structures(known_structs_dir)
    if not knowns:
        raise ValueError(f"No known-TPS structures found in {known_structs_dir}.")
    print(f"Designs: {len(designs)} ({dmode}); known refs: {len(knowns)} ({kmode})")

    panel = load_panel(panel_file) if panel_file else None
    if panel is not None:
        print(f"Loaded explicit SDR panel '{panel_file}' covering references: "
              f"{', '.join(panel.references())}")
    else:
        print(f"Using structure-derived SDR panel (residues within {panel_cutoff} A "
              "of the carboxylate-cage metal point of each neighbour).")

    nn = _rank1_neighbours(structural_topk, sequence_topk, set(knowns))
    print(f"Resolved rank-1 neighbours for {len(nn)}/{len(designs)} designs.")

    # Cache parsed neighbour structures (each known ref parsed at most once).
    neighbour_cache: Dict[str, Optional[ResidueInfo]] = {}

    def _neighbour_info(nid: str) -> Optional[ResidueInfo]:
        if nid in neighbour_cache:
            return neighbour_cache[nid]
        path = knowns.get(nid)
        info = None
        if path is not None:
            try:
                info = ResidueInfo(path)
            except Exception as exc:
                print(f"  [warn] could not parse neighbour '{nid}': {exc}")
        neighbour_cache[nid] = info
        return info

    rows: List[Dict[str, object]] = []
    n = len(designs)
    for i, (stem, path) in enumerate(designs.items(), start=1):
        row: Dict[str, object] = {"ID": stem}
        nn_entry = nn.get(stem)
        if nn_entry is None:
            row.update({
                "nearest_neighbour_id": "",
                "nearest_neighbour_similarity": np.nan,
                "similarity_space": "",
                "n_sdr_positions": 0,
                "sdr_identity": np.nan,
                "n_sdr_mismatches": np.nan,
                "specificity_divergence": False,
                "divergent_positions": "",
            })
        else:
            nid, sim, space = nn_entry
            if sim < tau_high:
                # Below the global floor: not the failure mode this metric targets.
                row.update({
                    "nearest_neighbour_id": nid,
                    "nearest_neighbour_similarity": float(sim),
                    "similarity_space": space,
                    "n_sdr_positions": 0,
                    "sdr_identity": np.nan,
                    "n_sdr_mismatches": np.nan,
                    "specificity_divergence": False,
                    "divergent_positions": "",
                })
            else:
                nbr_info = _neighbour_info(nid)
                try:
                    design_info = ResidueInfo(path)
                except Exception as exc:
                    print(f"  [warn] could not parse design '{stem}': {exc}")
                    design_info = None
                if nbr_info is None or design_info is None:
                    row.update({
                        "nearest_neighbour_id": nid,
                        "nearest_neighbour_similarity": float(sim),
                        "similarity_space": space,
                        "n_sdr_positions": 0,
                        "sdr_identity": np.nan,
                        "n_sdr_mismatches": np.nan,
                        "specificity_divergence": False,
                        "divergent_positions": "",
                    })
                else:
                    row.update(sdr_divergence_one(
                        design_info, nid, nbr_info, sim, space,
                        panel=panel, panel_cutoff=panel_cutoff,
                        map_tolerance=map_tolerance,
                        tau_high=tau_high, tau_low=tau_low,
                    ))
        rows.append(row)
        if i % 50 == 0 or i == n:
            print(f"  processed {i}/{n}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)
    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    n_flag = int(df["specificity_divergence"].sum())
    print(f"Wrote {len(df)} rows to {save_path} ({n_flag} flagged specificity_divergence).")
    return df
