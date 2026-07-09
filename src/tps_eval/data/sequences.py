from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple, Union

from Bio import SeqIO


SequenceRecord = Tuple[str, str]


def load_fasta_sequences(
    file_path: str,
    *,
    remove_padding: bool = True,
    load_identifiers: bool = False,
) -> Union[List[str], List[SequenceRecord]]:
    """Load sequences from a FASTA file.

    Args:
        file_path: Path to FASTA file.
        remove_padding: Remove '-' alignment gaps when True.
        load_identifiers: Return (id, sequence) tuples instead of sequences.
    """
    records = SeqIO.parse(file_path, "fasta")
    return read_sequences(
        records,
        remove_padding=remove_padding,
        load_identifiers=load_identifiers,
    )


def separate_identifiers(train: Sequence[SequenceRecord]) -> Tuple[List[str], List[str]]:
    identifiers = [record_id for record_id, _ in train]
    sequences = [sequence for _, sequence in train]
    return identifiers, sequences


def load_fasta_msa(file_path: str, *, load_identifiers: bool = False):
    return load_fasta_sequences(
        file_path,
        remove_padding=False,
        load_identifiers=load_identifiers,
    )


def read_sequences(
    records: Iterable,
    *,
    remove_padding: bool = True,
    load_identifiers: bool = False,
) -> Union[List[str], List[SequenceRecord]]:
    sequences: Union[List[str], List[SequenceRecord]] = []
    for record in records:
        seq = str(record.seq)
        if remove_padding:
            seq = seq.replace("-", "")
        if load_identifiers:
            sequences.append((str(record.id), seq))
        else:
            sequences.append(seq)
    return sequences
