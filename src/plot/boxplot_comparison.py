from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.load_results import load_results  # noqa: E402
from plot.constants import LOAD, MAX_VAL, MIN_VAL, THRESHOLD, TICKS  # noqa: E402
from plot.target_config import auto_ticks, resolve_range  # noqa: E402


def boxplot_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    target: str,
    *,
    dfs: Optional[Sequence] = None,
    save_dir: Optional[str] = None,
) -> None:
    """Boxplot of `target` across datasets.

    Expects each fasta file's auxiliary CSVs to live next to it. See
    `load_results` for details. Pass `dfs` (a per-dataset list of already-loaded
    DataFrames) to bypass fasta-based loading — used for structure targets,
    which are loaded from `<structs_dir>_<tool>.csv` files instead.
    """
    threshold = THRESHOLD.get(target)

    if dfs is None:
        dfs = [load_results(fp, load=LOAD[target]) for fp in fasta_paths]

    # When zero input CSVs survived (e.g. a structure-only run with no
    # sequence-branch outputs), the loaded frames carry no `target` column —
    # skip cleanly, mirroring the per-target "[skip] missing input" path.
    if not any(target in df.columns for df in dfs):
        print(f"  [skip] target {target}: no input data")
        return

    all_data: List[List[float]] = [
        [float(v) for v in df[target].dropna().tolist()]
        if target in df.columns
        else []
        for df in dfs
    ]

    # Fixed scale when defined in constants; otherwise auto-range from the data.
    min_val, max_val = resolve_range(target, MIN_VAL, MAX_VAL, all_data)
    ticks = TICKS.get(target)
    if ticks is None:
        ticks = auto_ticks(min_val, max_val)

    fig, ax = plt.subplots()

    positions = list(range(1, len(data_names) + 1))
    bp = ax.boxplot(
        all_data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        vert=True,
        medianprops={"color": "black"},
    )
    for patch, color in zip(bp["boxes"], data_colors):
        patch.set_facecolor(color)

    ax.set_title(target + " boxplot")
    ax.set_xticks(positions)
    ax.set_xticklabels(data_names)
    ax.set_ylabel(target)
    ax.set_yticks(ticks)
    ax.set_ylim(min_val, max_val)
    ax.yaxis.grid(True, linestyle="--", color="gray")

    if threshold is not None:
        ax.axhline(threshold, color="red", linestyle="--", linewidth=2)

    fig.tight_layout()

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, target + "_boxplot.png"))

    plt.close(fig)
