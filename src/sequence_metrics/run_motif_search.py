from __future__ import annotations

import sys

from motif_search import motif_search


def main() -> None:
    fasta_path = sys.argv[1]
    motifs = sys.argv[2:]
    motif_search(fasta_path, motifs, save=True)


if __name__ == "__main__":
    main()
