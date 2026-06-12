"""Codon optimization of a protein into a synthesis-ready coding sequence (CDS).

This is the FIRST step of order preparation. The protein is reverse-translated and
codon-optimized for a target organism (default: *Saccharomyces cerevisiae*). Three
correctness/synthesizability constraints are enforced:

1. **No internal BsaI/BsmBI sites** — those Type IIS sites are introduced deliberately and
   only by the overhangs (see ``overhangs.py``); a copy inside the CDS would fragment the
   part during Golden Gate assembly. (The previous web workflow let 2 sequences slip through
   with internal BsaI sites — this prevents that.)
2. **Homopolymer cap** — no single-nucleotide run longer than ``max_homopolymer`` (default 6);
   long runs cause synthesis slippage/indels and can act as cryptic yeast termination signals.
3. **GC window** — GC kept within ``[gc_min, gc_max]`` over every ``gc_window``-bp sliding
   window (local GC, which is what synthesis vendors check), avoiding extreme-GC stretches.

Constraints (1)-(3) are *hard*; codon usage is optimized as a soft objective within them.
The default method ``match_codon_usage`` samples codons to match the organism's natural
distribution (like the cad-sge web tool), which keeps GC and codon diversity natural;
``use_best_codon`` maximizes CAI but produces low-GC, repetitive, homopolymer-prone CDSs.

If the synthesizability constraints (2)-(3) make a particular protein infeasible (an
amino-acid stretch with only AT-rich synonymous codons can fight a GC floor), they are
**relaxed step by step** until a solution is found, and a note is appended to ``warnings``.
Translation and enzyme-site avoidance are NEVER relaxed (they are correctness-critical).

The stop codon is part of the optimized region (so the optimizer picks the organism's
preferred stop, e.g. TAA in yeast). Overhangs are added AFTERWARDS and never modified.

Backend: DNAChisel (https://edinburgh-genome-foundry.github.io/DnaChisel/) + the codon
tables from ``python_codon_tables``. Both are pip-installed in the ``tps_eval`` env
(see ``requirements.txt``).
"""
from __future__ import annotations

import random as _random

import numpy as _np
from dnachisel import (
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    EnzymeSitePattern,
    NoSolutionError,
    reverse_translate,
)

# ---------------------------------------------------------------------------
# Organism handling — friendly aliases -> python_codon_tables identifiers.
# Yeast is the only organism we need today; the map keeps the door open for others
# (any name/taxid accepted by python_codon_tables also works directly).
# ---------------------------------------------------------------------------
ORGANISM_ALIASES = {
    "yeast": "s_cerevisiae_4932",
    "s_cerevisiae": "s_cerevisiae_4932",
    "saccharomyces_cerevisiae": "s_cerevisiae_4932",
    "cerevisiae": "s_cerevisiae_4932",
    # examples for future use:
    "e_coli": "e_coli_316407",
    "escherichia_coli": "e_coli_316407",
    "h_sapiens": "h_sapiens_9606",
    "human": "h_sapiens_9606",
}
DEFAULT_ORGANISM = "yeast"

# Type IIS enzymes whose recognition sites must NOT appear inside the CDS. Both strands
# are handled automatically by DNAChisel's EnzymeSitePattern.
GOLDEN_GATE_ENZYMES = ("BsaI", "BsmBI")

# Codon-usage strategy and synthesizability defaults.
DEFAULT_METHOD = "match_codon_usage"   # sample to match natural codon distribution
DEFAULT_MAX_HOMOPOLYMER = 6            # longest allowed single-nucleotide run (0 disables)
DEFAULT_GC_MIN = 0.30                  # GC-window lower bound (fraction; None disables)
DEFAULT_GC_MAX = 0.65                  # GC-window upper bound (fraction)
DEFAULT_GC_WINDOW = 50                 # GC sliding-window size in bp (0 disables)
DEFAULT_SEED = 0                       # RNG seed for reproducible sampling (None = random)

_STOP_CODONS = ("TAA", "TAG", "TGA")


def resolve_organism(organism: str) -> str:
    """Map a friendly organism name to a python_codon_tables identifier (pass-through
    if it is already an identifier/taxid)."""
    return ORGANISM_ALIASES.get(organism.strip().lower(), organism)


def _build_constraints(seq, avoid_enzymes, max_homopolymer, gc_min, gc_max, gc_window):
    """Assemble the DNAChisel constraint list. Translation + enzyme-site avoidance are
    always present; the homopolymer cap and GC window are included only when enabled."""
    constraints = [EnforceTranslation()]
    constraints += [AvoidPattern(EnzymeSitePattern(enz)) for enz in avoid_enzymes]
    if max_homopolymer:
        # Forbidding a run of (cap + 1) of each base caps every homopolymer at `cap`.
        constraints += [AvoidPattern(f"{max_homopolymer + 1}x{base}") for base in "ATGC"]
    if gc_window and gc_min is not None and gc_max is not None:
        window = min(int(gc_window), len(seq))  # guard short CDSs
        constraints.append(EnforceGCContent(mini=gc_min, maxi=gc_max, window=window))
    return constraints


def _relaxation_ladder(max_homopolymer, gc_min, gc_max, gc_window):
    """Build the progressive relaxation ladder as (note, hp, gc_min, gc_max, gc_window)
    tuples, strictest first. Only loosens the synthesizability constraints, and only the
    steps that actually relax an enabled constraint."""
    ladder = [(None, max_homopolymer, gc_min, gc_max, gc_window)]
    gc_on = bool(gc_window) and gc_min is not None and gc_max is not None
    if gc_on:
        for delta in (0.05, 0.10):
            lo, hi = max(0.0, gc_min - delta), min(1.0, gc_max + delta)
            ladder.append((f"GC window widened to {lo:.2f}-{hi:.2f}",
                           max_homopolymer, lo, hi, gc_window))
        ladder.append(("GC-window constraint dropped", max_homopolymer, None, None, None))
    if max_homopolymer:
        suffix = " and GC-window constraint dropped" if gc_on else ""
        ladder.append((f"homopolymer cap raised to {max_homopolymer + 2}{suffix}",
                       max_homopolymer + 2, None, None, None))
        ladder.append(("homopolymer + GC-window constraints dropped (kept translation + "
                       "enzyme-site avoidance)", 0, None, None, None))
    return ladder


def codon_optimize(
    protein: str,
    organism: str = DEFAULT_ORGANISM,
    avoid_enzymes: tuple[str, ...] = GOLDEN_GATE_ENZYMES,
    add_stop: bool = True,
    method: str = DEFAULT_METHOD,
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER,
    gc_min: float | None = DEFAULT_GC_MIN,
    gc_max: float | None = DEFAULT_GC_MAX,
    gc_window: int = DEFAULT_GC_WINDOW,
    seed: int | None = DEFAULT_SEED,
    warnings: list[str] | None = None,
) -> str:
    """Reverse-translate and codon-optimize ``protein`` into a CDS for ``organism``.

    Parameters
    ----------
    protein:
        Amino-acid sequence (single-letter). A trailing ``*`` is tolerated and stripped.
    organism:
        Friendly name (see ``ORGANISM_ALIASES``) or a python_codon_tables identifier.
    avoid_enzymes:
        Type IIS enzymes whose sites are removed from the CDS (both strands). NEVER relaxed.
    add_stop:
        If True (default) a stop codon is appended and optimized with the CDS, so the
        optimizer chooses the organism-preferred stop (TAA in yeast).
    method:
        DNAChisel CodonOptimize method. ``match_codon_usage`` (default) matches the natural
        codon distribution; ``use_best_codon`` is deterministic but low-GC/repetitive.
    max_homopolymer:
        Longest allowed single-nucleotide run (default 6). 0 disables the cap.
    gc_min, gc_max:
        GC-window bounds as fractions (default 0.30-0.65). Set either to None to disable.
    gc_window:
        GC sliding-window size in bp (default 50). 0 disables.
    seed:
        RNG seed for reproducible sampling (default 0 → same input gives same output).
        Pass None for nondeterministic sampling. NOTE: seeds the global numpy/random RNGs.
    warnings:
        Optional list; any constraint-relaxation note is appended to it.

    Returns
    -------
    The optimized CDS (uppercase), starting with ``ATG`` and — if ``add_stop`` — ending in
    a stop codon, with no internal BsaI/BsmBI sites and (best effort) the homopolymer/GC
    constraints satisfied.
    """
    protein = protein.strip().upper().rstrip("*")
    if not protein:
        raise ValueError("Empty protein sequence.")

    if seed is not None:
        _np.random.seed(seed)
        _random.seed(seed)

    # Reverse-translate the protein, then append a placeholder stop the optimizer is free
    # to swap for the organism-preferred one (EnforceTranslation locks it to *a* stop).
    seq = reverse_translate(protein) + ("TAA" if add_stop else "")
    species = resolve_organism(organism)

    last_error: NoSolutionError | None = None
    for note, hp, gmin, gmax, gwin in _relaxation_ladder(max_homopolymer, gc_min, gc_max, gc_window):
        problem = DnaOptimizationProblem(
            sequence=seq,
            constraints=_build_constraints(seq, avoid_enzymes, hp, gmin, gmax, gwin),
            objectives=[CodonOptimize(species=species, method=method)],
            logger=None,
        )
        try:
            problem.resolve_constraints()
        except NoSolutionError as exc:
            last_error = exc
            continue
        problem.optimize()
        if note and warnings is not None:
            warnings.append(f"codon optimization relaxed: {note} (to find a feasible sequence)")
        return problem.sequence.upper()

    # Even the minimal set (translation + enzyme avoidance) failed — re-raise.
    raise last_error if last_error is not None else RuntimeError("codon optimization failed")
