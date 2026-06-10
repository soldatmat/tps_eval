from __future__ import annotations

import argparse

from motif_pair_distance import motif_pair_distance


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Residue (sequence) distance between the two class-I TPS "
        "metal-binding motifs (DDXXD-family and NSE/DTE) for every sequence in a "
        "FASTA. Writes a CSV keyed by ID next to the input "
        "(<input>_motif_pair_distance.csv); distances are NaN when a motif is absent."
    )
    parser.add_argument("fasta_path", help="FASTA file of sequences to evaluate.")
    args = parser.parse_args()
    motif_pair_distance(args.fasta_path, save=True)


if __name__ == "__main__":
    main()
