from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Iterable, List, Pattern, Sequence, Union

import pandas as pd

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers


def convert_to_regex(motif: str) -> Pattern[str]:
    return re.compile(motif)


def convert_motifs(motifs: Sequence[str]) -> List[Pattern[str]]:
    return [convert_to_regex(motif) for motif in motifs]


def load_df(fasta_path: str) -> pd.DataFrame:
    sequence_identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    return pd.DataFrame({"ID": sequence_identifiers, "sequence": sequences})


def find_motifs(df: pd.DataFrame, motifs: Iterable[Pattern[str]]) -> None:
    for motif in motifs:
        column_name = motif.pattern
        df[column_name] = [bool(motif.search(str(seq))) for seq in df["sequence"].tolist()]


def get_save_path(data_path: str) -> str:
    extension = data_path.split(".")[-1]
    return data_path[: -len(extension) - 1] + "_motifs.csv"


def motif_search(
    fasta_path: str,
    motifs: Union[Sequence[str], Sequence[Pattern[str]]],
    *,
    save: bool = True,
) -> pd.DataFrame:
    """Evaluate motif presence per sequence and optionally save a CSV report.

    Args:
        fasta_path: FASTA file containing evaluated sequences.
        motifs: Motif patterns as strings or compiled regex patterns.
        save: Save output CSV next to input FASTA when True.
    """
    if motifs and isinstance(motifs[0], str):
        motifs = convert_motifs(motifs)  # type: ignore[assignment]

    df = load_df(fasta_path)
    find_motifs(df, motifs)  # type: ignore[arg-type]

    if save:
        save_path = get_save_path(fasta_path)
        df.to_csv(save_path, index=False)

    return df
