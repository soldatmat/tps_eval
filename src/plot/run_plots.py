from __future__ import annotations

import sys
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from plot.plot_comparison import plot_comparison  # noqa: E402


def main(argv: list[str]) -> None:
    num_args = len(argv)
    if num_args < 3:
        raise SystemExit(
            "Invalid number of arguments. Expected at least 3: "
            "fasta_paths, data_names, data_colors."
        )

    fasta_paths = argv[0].split(",")
    data_names = argv[1].split(",")
    data_colors = argv[2].split(",")

    kwargs: dict = {}
    i = 3
    while i < num_args:
        flag = argv[i]
        if flag == "--targets":
            if i + 1 >= num_args:
                raise SystemExit("Missing value for --targets.")
            kwargs["targets"] = argv[i + 1].split(",")
            i += 2
        elif flag == "--save_dir":
            if i + 1 >= num_args:
                raise SystemExit("Missing value for --save_dir.")
            kwargs["save_dir"] = argv[i + 1]
            i += 2
        else:
            raise SystemExit(f"Unknown argument: {flag}")

    plot_comparison(fasta_paths, data_names, data_colors, **kwargs)


if __name__ == "__main__":
    main(sys.argv[1:])
