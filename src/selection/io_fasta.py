"""FASTA read/write helpers for the selection layer.

Reads reuse the repo's canonical loader (``data.sequences.load_fasta_sequences``); the
matching writer lives here because the selection layer is the first place that emits a
survivor subset FASTA (the metric tools only consume FASTAs, they never write them).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List

_SRC = Path(__file__).resolve().parent.parent  # the src/ dir
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from data.sequences import load_fasta_sequences  # noqa: E402


def read_fasta_map(fasta_path: str) -> Dict[str, str]:
    """Return {id: sequence} for a FASTA (padding stripped)."""
    records = load_fasta_sequences(fasta_path, load_identifiers=True)
    return {rid: seq for rid, seq in records}


def write_fasta(id_to_seq: Dict[str, str], output_path: str, order: List[str] = None,
                line_width: int = 0) -> None:
    """Write {id: seq} to a FASTA. ``order`` fixes record order (default: dict order);
    ``line_width`` > 0 wraps sequence lines, 0 writes each on one line."""
    ids = order if order is not None else list(id_to_seq)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as fh:
        for rid in ids:
            seq = id_to_seq[rid]
            fh.write(f">{rid}\n")
            if line_width and line_width > 0:
                for i in range(0, len(seq), line_width):
                    fh.write(seq[i:i + line_width] + "\n")
            else:
                fh.write(seq + "\n")
