from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Optional, Sequence

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from plot.boxplot_comparison import boxplot_comparison  # noqa: E402
from plot.constants import TARGETS  # noqa: E402
from plot.density_comparison import density_comparison  # noqa: E402


def plot_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    *,
    targets: Optional[Sequence[str]] = None,
    save_dir: Optional[str] = None,
) -> None:
    if targets is None:
        targets = TARGETS

    for target in targets:
        try:
            print(f"Generating boxplot for target: {target}...")
            boxplot_comparison(
                fasta_paths,
                data_names,
                data_colors,
                target,
                save_dir=save_dir,
            )

            print(f"Generating density plot for target: {target}...")
            density_comparison(
                fasta_paths,
                data_names,
                data_colors,
                target,
                direction="vertical",
                save_dir=save_dir,
            )
        except Exception as e:  # mirror Julia's per-target try/catch
            print(f"Error while plotting for target {target}: {e}")
            print("Stacktrace:")
            traceback.print_exc()
