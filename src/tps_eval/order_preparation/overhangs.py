"""Golden Gate / MoClo overhang library and CDS-wrapping logic.

This is the SECOND step of order preparation. The codon-optimized CDS is wrapped with
fixed flanking adapters so the part can be cloned by Golden Gate assembly. **The flanks
are constant and must never be codon-optimized or otherwise altered** — they carry the
BsaI/BsmBI recognition sites and fusion overhangs that the assembly relies on.

The library is transcribed from ``Plasmid_Generator.xlsx`` (the per-type ``F``/``J`` —
and for the Type 8 rows also ``G``/``I`` — columns). Every type shares the same constant
adapter scaffold:

    5' :  actcgacaac | CGTCTCa (BsmBI) | tcGGTCTCa (BsaI) | <left fusion>  ... CDS ...
    3' :  <right fusion> | tGAGACC (BsaI rev) | tGAGACG (BsmBI rev) | gttgtggtgt

The DEFAULT and only validated type is **Type 3** ("heterologous protein expression").
Its prefix/suffix were cross-checked against all 32 sequences of a previous real order
(``all_candidates_dna_fixed.txt``):

    prefix = actcgacaacCGTCTCatcGGTCTCaT
    suffix = ATCCtGAGACCtGAGACGgttgtggtgt   (preceded by the CDS's own stop codon)

Note from the sheet: *"Type 3 overhangs are only T on purpose. Don't touch them."*

The CDS-type rows bake a ``TAG`` stop into the right-fusion cell (Type 3 = ``TAGATCC``).
We strip that leading stop because the stop codon is produced by codon optimization
(and may be TAA/TGA), exactly as observed in the previous order.
"""
from __future__ import annotations

from dataclasses import dataclass

# Constant adapter scaffold shared by every overhang type (Plasmid_Generator.xlsx).
_FLANK5 = "actcgacaac"          # outer 5' flank (C column)
_BSMBI = "CGTCTCa"             # BsmBI (Esp3I) site + spacer (D column)
_BSAI = "tcGGTCTCa"           # BsaI site + spacer (E column)
_BSAI_R = "tGAGACC"            # BsaI site, reverse strand (K column)
_BSMBI_R = "tGAGACG"           # BsmBI site, reverse strand (L column)
_FLANK3 = "gttgtggtgt"         # outer 3' flank (M column)

_PREFIX_BASE = _FLANK5 + _BSMBI + _BSAI       # everything 5' of the left fusion overhang
_SUFFIX_BASE = _BSAI_R + _BSMBI_R + _FLANK3   # everything 3' of the right fusion overhang

_STOP_CODONS = ("TAA", "TAG", "TGA")


@dataclass(frozen=True)
class OverhangType:
    """Variable fusion overhangs of one MoClo/YTK part type (the constant scaffold is
    shared). ``left_fusion`` sits just 5' of the CDS; ``right_fusion`` just 3' of it
    (after the CDS stop codon, with any leading stop stripped at assembly time)."""

    name: str
    left_fusion: str
    right_fusion: str
    description: str = ""


# Transcribed from Plasmid_Generator.xlsx. Only Type 3 is validated against a real order.
OVERHANGS: dict[str, OverhangType] = {
    "Type 3": OverhangType("Type 3", "T", "TAGATCC", "heterologous protein expression (DEFAULT)"),
    "Type 1": OverhangType("Type 1", "CCCT", "AACG"),
    "Type 2": OverhangType("Type 2", "AACG", "TATG"),
    "Type 3a": OverhangType("Type 3a", "T", "TTCT", "CDS with C-terminal fusion (no stop)"),
    "Type 3b": OverhangType("Type 3b", "TTCT", "TAGATCC", "CDS with N-terminal fusion"),
    "Type 4": OverhangType("Type 4", "ATCC", "GCTG"),
    "Type 4a": OverhangType("Type 4a", "ATCC", "TGGC"),
    "Type 4b": OverhangType("Type 4b", "TGGC", "GCTG"),
    "Type 5": OverhangType("Type 5", "GCTG", "TACA"),
    "Type 6 plasmid": OverhangType("Type 6 plasmid", "TACA", "GAGT"),
    "Type 7 plasmid origin": OverhangType("Type 7 plasmid origin", "GAGT", "CCGA"),
    "Type 8 plasmid bacteria": OverhangType("Type 8 plasmid bacteria", "CCGAgcggccgc", "gcggccgcCCCT"),
    "Type 6 integration": OverhangType("Type 6 integration", "TACA", "GAGT"),
    "Type 7 integration": OverhangType("Type 7 integration", "GAGT", "CCGA"),
    "Type 8a integration": OverhangType("Type 8a integration", "CCGAgcggccgc", "gcggccgcCAAT"),
    "Type 8b integration": OverhangType("Type 8b integration", "CAAT", "CCCT"),
}

DEFAULT_OVERHANG = "Type 3"


def _strip_leading_stop(fusion: str) -> str:
    """Drop a leading in-frame stop codon from a right-fusion (the stop comes from the
    codon-optimized CDS, not the overhang)."""
    if fusion[:3].upper() in _STOP_CODONS:
        return fusion[3:]
    return fusion


def get_overhangs(overhang_type: str = DEFAULT_OVERHANG) -> tuple[str, str]:
    """Return ``(prefix, suffix)`` for ``overhang_type``.

    ``prefix`` is prepended to the CDS; ``suffix`` is appended directly after the CDS's
    stop codon. For Type 3 this yields the previous order's exact flanks:
    ``actcgacaacCGTCTCatcGGTCTCaT`` and ``ATCCtGAGACCtGAGACGgttgtggtgt``.
    """
    try:
        oh = OVERHANGS[overhang_type]
    except KeyError:
        raise KeyError(
            f"Unknown overhang type {overhang_type!r}. "
            f"Available: {', '.join(OVERHANGS)}"
        )
    prefix = _PREFIX_BASE + oh.left_fusion
    suffix = _strip_leading_stop(oh.right_fusion) + _SUFFIX_BASE
    return prefix, suffix


def add_overhangs(cds: str, overhang_type: str = DEFAULT_OVERHANG) -> str:
    """Wrap a codon-optimized CDS with the fixed flanks of ``overhang_type``.

    The CDS is expected to start with ``ATG`` and end in a stop codon. The flanks are
    constant and are never modified.
    """
    prefix, suffix = get_overhangs(overhang_type)
    return prefix + cds + suffix
