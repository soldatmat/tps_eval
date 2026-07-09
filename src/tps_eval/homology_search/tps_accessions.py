from __future__ import annotations

"""Shared TPS-vs-non-TPS classification core for the broad-homology-search tools.

Both the sequence search (swissprot_search) and the structure search
(foldseek_swissprot_search) map each database hit to a UniProt accession, then
classify the hit as a terpene synthase (TPS) or not by membership in a precomputed
set of Swiss-Prot TPS accessions.

The accession set lives in a committable data file (one accession per line),
generated ONCE from the UniProt REST API with the query::

    (reviewed:true) AND ((ec:4.2.3.*) OR (ec:5.5.1.*))

i.e. reviewed (Swiss-Prot) entries that are class I terpene synthases/cyclases
(EC 4.2.3.*) or class II terpene synthases/cyclases (EC 5.5.1.*).

IMPORTANT classification nuance: prenyltransferases / isoprenyl-diphosphate
synthases (EC 2.5.1.*) are deliberately NOT included — they are exactly the
related-but-different enzymes we want to be able to flag when a design's closest
relative turns out to be a non-TPS. So a hit to e.g. farnesyl-diphosphate synthase
correctly classifies as non-TPS.
"""

import os
from functools import lru_cache
from typing import Set

# UniProt accessions are alphanumeric; isoform suffixes ("-2") and version
# suffixes are stripped to the canonical accession before lookup.


def _canonical_accession(accession: str) -> str:
    """Strip isoform/version suffixes to the canonical UniProt accession.

    e.g. "P12345-2" -> "P12345", "P12345.1" -> "P12345".
    """
    acc = accession.strip()
    # Isoform suffix.
    acc = acc.split("-", 1)[0]
    # Version suffix (rare in our inputs but harmless).
    acc = acc.split(".", 1)[0]
    return acc


@lru_cache(maxsize=None)
def load_tps_accessions(path: str) -> frozenset:
    """Load the committable TPS accession file into a frozenset (one acc/line).

    Cached per path so both tools (and repeated calls) share a single read.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"TPS accession file not found: {path}. Point TPS_ACCESSIONS in paths.sh "
            "at src/homology_search/tps_uniprot_accessions.txt (regenerate via the UniProt "
            "REST query documented in src/homology_search/tps_accessions.py)."
        )
    accs: Set[str] = set()
    with open(path) as fh:
        for line in fh:
            acc = _canonical_accession(line)
            if acc:
                accs.add(acc)
    return frozenset(accs)


def is_tps(accession: str, tps_accessions: frozenset) -> bool:
    """True iff the (canonicalized) accession is in the TPS set."""
    return _canonical_accession(accession) in tps_accessions
