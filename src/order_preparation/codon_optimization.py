"""Codon optimization of a protein into a synthesis-ready coding sequence (CDS).

This is the FIRST step of order preparation. The protein is reverse-translated and
codon-optimized for a target organism (default: *Saccharomyces cerevisiae*). Crucially,
the recognition sites of the Golden Gate / MoClo Type IIS enzymes (BsaI, BsmBI) are
*removed from the CDS* — those sites are introduced deliberately and only by the
overhangs (see ``overhangs.py``); any copy inside the CDS would fragment the part during
assembly. The previous order had 2 sequences with internal BsaI sites that slipped
through the manual (web) workflow — this step prevents that class of defect.

The stop codon is part of the optimized region (so the optimizer picks the organism's
preferred stop, e.g. TAA in yeast). Overhangs are added AFTERWARDS and never modified.

Backend: DNAChisel (https://edinburgh-genome-foundry.github.io/DnaChisel/) + the codon
tables from ``python_codon_tables``. Both are pip-installed in the ``tps_eval`` env
(see ``requirements.txt``).
"""
from __future__ import annotations

from dnachisel import (
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
    EnzymeSitePattern,
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

_STOP_CODONS = ("TAA", "TAG", "TGA")


def resolve_organism(organism: str) -> str:
    """Map a friendly organism name to a python_codon_tables identifier (pass-through
    if it is already an identifier/taxid)."""
    return ORGANISM_ALIASES.get(organism.strip().lower(), organism)


def codon_optimize(
    protein: str,
    organism: str = DEFAULT_ORGANISM,
    avoid_enzymes: tuple[str, ...] = GOLDEN_GATE_ENZYMES,
    add_stop: bool = True,
    method: str = "use_best_codon",
) -> str:
    """Reverse-translate and codon-optimize ``protein`` into a CDS for ``organism``.

    Parameters
    ----------
    protein:
        Amino-acid sequence (single-letter). A trailing ``*`` is tolerated and stripped.
    organism:
        Friendly name (see ``ORGANISM_ALIASES``) or a python_codon_tables identifier.
    avoid_enzymes:
        Type IIS enzymes whose sites are removed from the CDS (both strands).
    add_stop:
        If True (default) a stop codon is appended and optimized together with the CDS,
        so the optimizer chooses the organism-preferred stop (TAA in yeast).
    method:
        DNAChisel CodonOptimize method. ``use_best_codon`` (default) is deterministic.

    Returns
    -------
    The optimized CDS (uppercase), starting with ``ATG`` and — if ``add_stop`` — ending
    in a stop codon, with no internal BsaI/BsmBI sites.
    """
    protein = protein.strip().upper().rstrip("*")
    if not protein:
        raise ValueError("Empty protein sequence.")

    # Reverse-translate the protein, then append a placeholder stop the optimizer is
    # free to swap for the organism-preferred one (EnforceTranslation locks it to *a*
    # stop; CodonOptimize picks the best/most-frequent stop codon).
    seq = reverse_translate(protein) + ("TAA" if add_stop else "")

    constraints = [EnforceTranslation()]
    constraints += [AvoidPattern(EnzymeSitePattern(enz)) for enz in avoid_enzymes]

    problem = DnaOptimizationProblem(
        sequence=seq,
        constraints=constraints,
        objectives=[CodonOptimize(species=resolve_organism(organism), method=method)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence.upper()
