from __future__ import annotations

"""Shared localization of the two class-I terpene-synthase metal-binding motifs.

Class-I TPS bind their Mg2+/Mn2+ cluster via two motifs on opposite walls of the
active-site cleft:

* the aspartate-rich **DDXXD** motif (here matched by the relaxed acidic family
  ``[DE][DE]..[DE]`` so the conservative D->E substitutions DExxD / EDxxD / EExxE
  still count — generative models frequently emit those, and they still supply the
  carboxylates that coordinate the first metal pair); and
* the **NSE/DTE** motif (the second metal-binding triad), consensus roughly
  ``(N/D)Dxx(S/T)xxxE``. We match it with ``(N|D)D(L|I|V).(S|T)...E`` — i.e. the
  literature core ``(N/D)D``...``(S/T)``...``E`` with the common hydrophobic
  ``[LIV]`` at position 3 — the SAME regex already used by ``run_motif_search``'s
  default motif list, so the boolean motif_search columns and these positions agree.

This module returns, for a sequence, WHERE each motif's best match sits (so the
sequence and structural distance tools share one source of truth). Position
indices are reported as **0-based** ``start``/``end`` half-open spans
(``seq[start:end]`` is the matched substring) — convenient for slicing and for
indexing into a Biopython residue list — plus a 1-based ``start_1`` for human/CSV
readability. "Best" match = the FIRST match scanning left-to-right (motifs are
short and a well-formed class-I TPS has one of each in a canonical order; when a
sequence has several, the first is the conventional pick and keeps the two tools
deterministic and mutually consistent).
"""

import re
from dataclasses import dataclass
from typing import List, Optional, Pattern

# Relaxed acidic DDXXD family (strict DD..D is a subset). Matching the widest of
# the three nested variants from run_motif_search so we localize the motif even
# when a model emitted a D->E-substituted first metal site.
DDXXD_PATTERN: Pattern[str] = re.compile(r"[DE][DE]..[DE]")

# NSE/DTE second metal-binding motif. Identical to the run_motif_search default.
NSE_DTE_PATTERN: Pattern[str] = re.compile(r"(N|D)D(L|I|V).(S|T)...E")

# Within each motif, the residues whose side chains actually coordinate the metals
# (used by the structural tool to pick which CA atoms span the metal cluster):
#   DDXXD family: the three acidic positions are indices 0, 1 and 4 of the 5-mer.
#   NSE/DTE     : the (N/D), the following D, the (S/T) and the terminal E, i.e.
#                 indices 0, 1, 4 and 8 of the 9-mer.
DDXXD_COORDINATING_OFFSETS = (0, 1, 4)
NSE_DTE_COORDINATING_OFFSETS = (0, 1, 4, 8)


@dataclass
class MotifMatch:
    """A localized motif hit. ``start``/``end`` are 0-based half-open
    (``sequence[start:end] == matched``); ``start_1`` is 1-based."""

    matched: str
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive (== start + len(matched))
    start_1: int        # 1-based start, for CSV/human use

    @property
    def coordinating_offsets(self) -> tuple:  # overridden per-motif at build time
        return tuple(range(len(self.matched)))


def _first_match(sequence: str, pattern: Pattern[str]) -> Optional[re.Match]:
    return pattern.search(sequence) if sequence else None


def locate_ddxxd(sequence: str) -> Optional[MotifMatch]:
    """First DDXXD-family (``[DE][DE]..[DE]``) match, or None."""
    m = _first_match(sequence, DDXXD_PATTERN)
    if m is None:
        return None
    return MotifMatch(m.group(0), m.start(), m.end(), m.start() + 1)


def locate_nse_dte(sequence: str) -> Optional[MotifMatch]:
    """First NSE/DTE (``(N|D)D(L|I|V).(S|T)...E``) match, or None."""
    m = _first_match(sequence, NSE_DTE_PATTERN)
    if m is None:
        return None
    return MotifMatch(m.group(0), m.start(), m.end(), m.start() + 1)


def coordinating_indices(match: MotifMatch, offsets: tuple) -> List[int]:
    """0-based sequence indices of the metal-coordinating residues of ``match``,
    given the motif's coordinating-residue offsets. Offsets that fall outside the
    matched span (should not happen for the fixed-length patterns) are dropped."""
    return [match.start + off for off in offsets if 0 <= off < len(match.matched)]
