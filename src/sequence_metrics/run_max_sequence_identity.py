from __future__ import annotations

import sys

from max_sequence_identity import max_sequence_identity


def main() -> None:
    num_args = len(sys.argv) - 1
    if num_args == 1:
        fasta_path = sys.argv[1]
        max_sequence_identity(fasta_path)
    elif num_args == 2:
        fasta_path = sys.argv[1]
        train_path = sys.argv[2]
        max_sequence_identity(fasta_path, train_path=train_path)
    else:
        raise ValueError(f"Invalid number of arguments. Expected 1 or 2, got {num_args}.")


if __name__ == "__main__":
    main()
