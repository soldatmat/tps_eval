"""Order preparation: amino-acid designs -> codon-optimized, overhang-flanked DNA.

Pipeline (per design):
    0. normalize the amino-acid termini so the design is a complete ORF — prepend a start
       Met (-> ATG start codon) and ensure a C-terminal stop codon are present, ADDING
       either if missing                                 (``normalize_termini``)
    1. codon-optimize the protein into a CDS for the target organism, removing internal
       BsaI/BsmBI sites          (``codon_optimization.codon_optimize``)
    2. wrap the CDS with the fixed Golden Gate flanks  (``overhangs.add_overhangs``)
    3. validate the assembled construct               (``validate_construct``)

This is NOT part of the evaluation / submit-all pipeline — it is a standalone ordering
utility. It runs on a login node (fast, no SLURM).

Input  : a FASTA of protein sequences, or a CSV with id + amino-acid columns
         (auto-detected by extension), or a single inline sequence.
Output : ``<prefix>_order.csv``  — id, status, start_added, stop_added, protein, cds,
                                   ordered_sequence, length_nt, warnings
         ``<prefix>_order.txt``  — ``id,ordered_sequence`` per line (matches the format of
                                   the previous order's all_candidates_dna_fixed.txt)
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from tps_eval.order_preparation.codon_optimization import (
    DEFAULT_GC_MAX,
    DEFAULT_GC_MIN,
    DEFAULT_GC_WINDOW,
    DEFAULT_MAX_HOMOPOLYMER,
    DEFAULT_METHOD,
    DEFAULT_ORGANISM,
    DEFAULT_SEED,
    GOLDEN_GATE_ENZYMES,
    codon_optimize,
)
from tps_eval.order_preparation.overhangs import DEFAULT_OVERHANG, add_overhangs, get_overhangs

_FASTA_EXTS = {".fasta", ".fa", ".faa", ".fas", ".fna"}
_TABLE_EXTS = {".csv", ".tsv"}
_STOP_CODONS = ("TAA", "TAG", "TGA")

# Recognition sequences (5'->3') of the Golden Gate enzymes, for validation.
_ENZYME_SITES = {"BsaI": "GGTCTC", "BsmBI": "CGTCTC"}

# Re-optimization attempts to clear an unintended Type IIS site before failing a design.
DEFAULT_TYPEIIS_RETRIES = 5

_ID_COLUMNS = ("id", "name", "identifier", "seq_id", "sequence_id", "design_id")
_SEQ_COLUMNS = (
    "sequence", "aa", "aa_sequence", "amino_acid_sequence", "aa_seq",
    "protein", "protein_sequence", "seq",
)


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _max_homopolymer_run(s: str) -> int:
    """Length of the longest single-nucleotide run in ``s`` (0 for empty)."""
    return max((len(m.group(0)) for m in re.finditer(r"(.)\1*", s)), default=0)


def _gc_window_extremes(s: str, window: int) -> tuple[float, float]:
    """Min and max GC fraction over every ``window``-bp sliding window of ``s``. If ``s``
    is shorter than the window, returns its overall GC for both."""
    s = s.upper()
    n = len(s)
    prefix = [0] * (n + 1)
    for i, ch in enumerate(s):
        prefix[i + 1] = prefix[i] + (1 if ch in "GC" else 0)
    if n <= window:
        gc = prefix[n] / n if n else 0.0
        return gc, gc
    lo, hi = 1.0, 0.0
    for i in range(0, n - window + 1):
        gc = (prefix[i + window] - prefix[i]) / window
        lo, hi = min(lo, gc), max(hi, gc)
    return lo, hi


def _typeiis_violations(
    full: str, prefix: str, suffix: str, avoid_enzymes: tuple[str, ...]
) -> list[str]:
    """Type IIS recognition sites that are NOT fully inside a flank — i.e. inside the CDS
    or spanning a CDS<->flank junction. Empty list == only the deliberate flank sites are
    present. These are FATAL for Golden Gate assembly (the part would be cut)."""
    violations: list[str] = []
    cds_start, cds_end = len(prefix), len(full) - len(suffix)
    up = full.upper()
    for enz in avoid_enzymes:
        site = _ENZYME_SITES.get(enz)
        if site is None:
            continue
        for pattern in (site, _revcomp(site)):
            for m in re.finditer(pattern, up):
                s, e = m.start(), m.end()
                if not (e <= cds_start or s >= cds_end):  # not wholly inside a flank
                    violations.append(f"{enz} site at nt {s} overlaps the CDS/junction")
    return violations


# ---------------------------------------------------------------------------
# Termini normalization (complete-ORF guarantee)
# ---------------------------------------------------------------------------
# The amino acid whose codon is the start codon (ATG). A design must begin with it for its
# CDS to open with a start codon.
_START_AA = "M"


def normalize_termini(protein: str) -> tuple[str, bool, bool]:
    """Ensure the amino-acid design encodes a complete ORF: an N-terminal start (Met, i.e.
    an ``ATG`` start codon) and a C-terminal stop codon. Normalizes both *at the amino-acid
    level*, before codon optimization, so each fix lands in the correct place — the start
    codon at the very start of the CDS (just after the 5' overhang) and the stop codon at
    the very end of the CDS (just before the 3' overhang).

    * **Start** — if the sequence does not begin with ``M``, an ``M`` is prepended; its
      codon becomes the ``ATG`` start codon.
    * **Stop** — a trailing ``*`` (if present) marks a design-supplied stop and is stripped
      here; the stop codon itself is emitted by codon optimization (``add_stop``), which
      lets the optimizer pick the organism-preferred stop. A design without a trailing
      ``*`` is treated as lacking a stop, so one is being added.

    Returns ``(protein, start_added, stop_added)``: the normalized amino-acid sequence
    (Met-prefixed as needed, with any trailing ``*`` removed) plus booleans recording
    whether a start / stop codon had to be added for this design.
    """
    protein = protein.strip().upper()
    had_stop = protein.endswith("*")
    protein = protein.rstrip("*")
    if not protein:
        raise ValueError("Empty protein sequence.")
    start_added = not protein.startswith(_START_AA)
    if start_added:
        protein = _START_AA + protein
    stop_added = not had_stop
    return protein, start_added, stop_added


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------
def load_designs(
    input_path: str | Path,
    id_column: str | None = None,
    seq_column: str | None = None,
) -> list[tuple[str, str]]:
    """Load ``(id, protein)`` pairs from a FASTA or CSV/TSV file (auto-detected)."""
    path = Path(input_path)
    ext = path.suffix.lower()
    if ext in _FASTA_EXTS:
        return _load_fasta(path)
    if ext in _TABLE_EXTS:
        return _load_table(path, id_column, seq_column)
    # Unknown extension: try FASTA, fall back to table.
    try:
        designs = _load_fasta(path)
        if designs:
            return designs
    except Exception:
        pass
    return _load_table(path, id_column, seq_column)


def _load_fasta(path: Path) -> list[tuple[str, str]]:
    return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(str(path), "fasta")]


def _load_table(path: Path, id_column: str | None, seq_column: str | None) -> list[tuple[str, str]]:
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(path, sep=sep)
    lower = {c.lower(): c for c in df.columns}
    id_col = id_column or next((lower[c] for c in _ID_COLUMNS if c in lower), None)
    seq_col = seq_column or next((lower[c] for c in _SEQ_COLUMNS if c in lower), None)
    if seq_col is None:
        raise ValueError(
            f"Could not find an amino-acid column in {path.name} (columns: "
            f"{list(df.columns)}). Pass --seq-column explicitly."
        )
    if id_col is None:
        # No id column -> use the row index.
        return [(str(i), str(s)) for i, s in enumerate(df[seq_col])]
    return [(str(i), str(s)) for i, s in zip(df[id_col], df[seq_col])]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_construct(
    full: str,
    protein: str,
    overhang_type: str = DEFAULT_OVERHANG,
    avoid_enzymes: tuple[str, ...] = GOLDEN_GATE_ENZYMES,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    gc_min: float | None = DEFAULT_GC_MIN,
    gc_max: float | None = DEFAULT_GC_MAX,
    gc_window: int = DEFAULT_GC_WINDOW,
) -> list[str]:
    """Independently re-check an assembled construct. Returns a list of human-readable
    warnings (empty == all checks passed).

    Checks: starts with the expected prefix; CDS starts with ATG; CDS length is a
    multiple of 3; CDS translates back to ``protein`` and ends in a stop codon; ends with
    the expected suffix; no BsaI/BsmBI site lies inside the CDS or across a junction (only
    the deliberate sites in the flanks are allowed); no homopolymer run in the full
    construct exceeds ``max_homopolymer``; and the CDS GC stays within ``[gc_min, gc_max]``
    over every ``gc_window``-bp window. The homopolymer/GC checks mirror the optimization
    targets, so a relaxed sequence will (intentionally) flag here too.
    """
    warnings: list[str] = []
    prefix, suffix = get_overhangs(overhang_type)
    protein = protein.strip().upper().rstrip("*")

    if not full.startswith(prefix):
        warnings.append(f"does not start with the {overhang_type} prefix")
    if not full.endswith(suffix):
        warnings.append(f"does not end with the {overhang_type} suffix")

    cds = full[len(prefix): len(full) - len(suffix)]
    if not cds.upper().startswith("ATG"):
        warnings.append(f"CDS does not start with a start codon (ATG): {cds[:3]!r}")
    if len(cds) % 3 != 0:
        warnings.append(f"CDS length {len(cds)} is not a multiple of 3")
    if cds[-3:].upper() not in _STOP_CODONS:
        warnings.append(f"CDS does not end with a stop codon: {cds[-3:]!r}")
    else:
        translated = str(Seq(cds).translate(to_stop=False)).rstrip("*")
        if translated != protein:
            warnings.append("translated CDS does not match the input protein")

    # No Type IIS site may overlap the CDS region (allowed only fully inside a flank).
    warnings += _typeiis_violations(full, prefix, suffix, avoid_enzymes)

    # Homopolymer cap — checked on the full construct (junctions included).
    up = full.upper()
    if max_homopolymer:
        run = _max_homopolymer_run(up)
        if run > max_homopolymer:
            warnings.append(f"homopolymer run of {run} nt exceeds cap of {max_homopolymer}")

    # GC window — checked on the CDS (the region we control; flanks are fixed).
    if gc_window and gc_min is not None and gc_max is not None:
        lo, hi = _gc_window_extremes(cds, min(int(gc_window), len(cds)))
        if lo < gc_min - 1e-9 or hi > gc_max + 1e-9:
            warnings.append(
                f"CDS GC over {gc_window}-bp windows is {lo*100:.0f}-{hi*100:.0f}%, "
                f"outside [{gc_min*100:.0f}, {gc_max*100:.0f}]%"
            )
    return warnings


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def prepare_one(
    protein: str,
    organism: str = DEFAULT_ORGANISM,
    overhang_type: str = DEFAULT_OVERHANG,
    avoid_enzymes: tuple[str, ...] = GOLDEN_GATE_ENZYMES,
    method: str = DEFAULT_METHOD,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    gc_min: float | None = DEFAULT_GC_MIN,
    gc_max: float | None = DEFAULT_GC_MAX,
    gc_window: int = DEFAULT_GC_WINDOW,
    seed: int | None = DEFAULT_SEED,
    max_attempts: int = DEFAULT_TYPEIIS_RETRIES,
) -> dict:
    """Run the full pipeline for a single protein and return a result row.

    An unintended Type IIS (BsaI/BsmBI) site inside the CDS or across a junction is FATAL
    for Golden Gate assembly. Because a junction site depends on the sampled codons, we
    re-optimize with a fresh seed up to ``max_attempts`` times to clear it; if it still
    persists the design is returned with ``status="FAILED"`` (the caller excludes it from
    the order .txt and exits non-zero). Homopolymer/GC issues are SOFT — they remain
    warnings and never fail a design.

    Before optimizing, the termini are normalized to a complete ORF (``normalize_termini``):
    a start Met (-> ATG) is prepended if missing and a stop codon is ensured, and the
    returned dict records whether each was ADDED (``start_added`` / ``stop_added``) — with
    a matching warning so a silent addition never happens.

    Returns a dict with ``status`` ("ok" | "FAILED"), ``start_added``, ``stop_added``,
    ``protein``, ``cds``, ``ordered_sequence``, ``length_nt``, ``warnings``.
    """
    prefix, suffix = get_overhangs(overhang_type)
    protein_norm, start_added, stop_added = normalize_termini(protein)
    termini_warnings: list[str] = []
    if start_added:
        termini_warnings.append("start codon (ATG) added: design did not begin with Met (M)")
    if stop_added:
        termini_warnings.append("stop codon added at the C-terminus: design carried no terminal stop")
    attempts = max(1, max_attempts)
    full = cds = ""
    violations: list[str] = []
    for attempt in range(attempts):
        opt_warnings: list[str] = []  # codon-optimization relaxation notes
        # Vary the seed per attempt so a fresh draw clears a junction site (None stays
        # nondeterministic, advancing the global RNG each call).
        attempt_seed = None if seed is None else seed + attempt
        cds = codon_optimize(
            protein_norm, organism=organism, avoid_enzymes=avoid_enzymes, method=method,
            max_homopolymer=max_homopolymer, gc_min=gc_min, gc_max=gc_max,
            gc_window=gc_window, seed=attempt_seed, warnings=opt_warnings,
        )
        full = add_overhangs(cds, overhang_type)
        violations = _typeiis_violations(full, prefix, suffix, avoid_enzymes)
        if violations:
            continue  # fatal site present — re-optimize with a new seed

        warnings = termini_warnings + opt_warnings + validate_construct(
            full, protein_norm, overhang_type, avoid_enzymes,
            max_homopolymer=max_homopolymer, gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
        )
        if attempt:
            warnings.insert(0, f"cleared a Type IIS site by re-optimizing (attempt {attempt + 1})")
        return {
            "status": "ok",
            "start_added": start_added,
            "stop_added": stop_added,
            "protein": protein_norm,
            "cds": cds,
            "ordered_sequence": full,
            "length_nt": len(full),
            "warnings": "; ".join(warnings),
        }

    # Exhausted all attempts — keep the last sequence for inspection but mark it FAILED.
    fatal = (f"FATAL: unresolved Type IIS site after {attempts} re-optimization "
             f"attempt(s) ({'; '.join(violations)})")
    return {
        "status": "FAILED",
        "start_added": start_added,
        "stop_added": stop_added,
        "protein": protein_norm,
        "cds": cds,
        "ordered_sequence": full,
        "length_nt": len(full),
        "warnings": "; ".join(termini_warnings + [fatal]),
    }


def prepare_order(
    input_path: str | Path,
    output_prefix: str | Path | None = None,
    organism: str = DEFAULT_ORGANISM,
    overhang_type: str = DEFAULT_OVERHANG,
    id_column: str | None = None,
    seq_column: str | None = None,
    method: str = DEFAULT_METHOD,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    gc_min: float | None = DEFAULT_GC_MIN,
    gc_max: float | None = DEFAULT_GC_MAX,
    gc_window: int = DEFAULT_GC_WINDOW,
    seed: int | None = DEFAULT_SEED,
    max_attempts: int = DEFAULT_TYPEIIS_RETRIES,
    save: bool = True,
) -> pd.DataFrame:
    """Prepare an ordering table from a file of protein designs.

    Writes ``<output_prefix>_order.csv`` (all designs, with a ``status`` column) and
    ``<output_prefix>_order.txt`` (only the ``status == "ok"`` designs) when ``save`` is
    True. Designs with an unresolved Type IIS site are marked ``FAILED`` and EXCLUDED from
    the .txt. ``output_prefix`` defaults to the input path without its suffix. Returns the
    result DataFrame (inspect ``df["status"]`` to detect failures).
    """
    designs = load_designs(input_path, id_column, seq_column)
    if not designs:
        raise ValueError(f"No sequences found in {input_path}.")

    rows = []
    for design_id, protein in designs:
        row = prepare_one(
            protein, organism, overhang_type, method=method,
            max_homopolymer=max_homopolymer, gc_min=gc_min, gc_max=gc_max,
            gc_window=gc_window, seed=seed, max_attempts=max_attempts,
        )
        row = {"id": design_id, **row}
        rows.append(row)
        if row["status"] != "ok":
            print(f"  [FAIL] {design_id}: {row['warnings']}")
        elif row["warnings"]:
            print(f"  [WARN] {design_id}: {row['warnings']}")
    df = pd.DataFrame(
        rows,
        columns=["id", "status", "start_added", "stop_added", "protein", "cds",
                 "ordered_sequence", "length_nt", "warnings"],
    )

    # Report termini additions up front (how many out of how many for each), so a silent
    # start/stop insertion never happens.
    n = len(df)
    n_start_added = int(df["start_added"].sum())
    n_stop_added = int(df["stop_added"].sum())
    if n_start_added:
        print(f"  [WARN] added a start codon (ATG) to {n_start_added}/{n} design(s) "
              f"that did not begin with Met (M).")
    if n_stop_added:
        print(f"  [WARN] added a stop codon to {n_stop_added}/{n} design(s) "
              f"that carried no terminal stop.")

    n_failed = int((df["status"] != "ok").sum())
    n_warn = int(((df["status"] == "ok") & (df["warnings"] != "")).sum())
    print(
        f"Prepared {len(df)} construct(s) for organism={organism!r}, overhang={overhang_type!r}. "
        f"{n_failed} FAILED, {n_warn} with warnings."
    )

    if save:
        prefix = Path(output_prefix) if output_prefix else Path(input_path).with_suffix("")
        csv_path = prefix.with_name(prefix.name + "_order.csv")
        txt_path = prefix.with_name(prefix.name + "_order.txt")
        df.to_csv(csv_path, index=False)  # all designs, incl. FAILED, with status
        ok = df[df["status"] == "ok"]
        with open(txt_path, "w") as fh:
            for _, r in ok.iterrows():
                fh.write(f"{r['id']},{r['ordered_sequence']}\n")
        print(f"Wrote {csv_path} ({len(df)} rows incl. status)")
        print(f"Wrote {txt_path} ({len(ok)} order-ready sequence(s))")
        if n_failed:
            print(
                f"  !! {n_failed} design(s) FAILED (unresolved Type IIS site) and were "
                f"EXCLUDED from {txt_path.name}; see the status/warnings columns in {csv_path.name}."
            )
    return df
