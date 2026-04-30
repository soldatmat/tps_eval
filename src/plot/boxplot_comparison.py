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


def boxplot_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    target: str,
    *,
    save_dir: Optional[str] = None,
) -> None:
    """Boxplot of `target` across datasets.

    Expects each fasta file's auxiliary CSVs to live next to it. See
    `load_results` for details.
    """
    load_list = LOAD[target]
    min_val = MIN_VAL[target]
    max_val = MAX_VAL[target]
    ticks = TICKS[target]
    threshold = THRESHOLD[target]

    all_dfs = [load_results(fp, load=load_list) for fp in fasta_paths]
    all_data: List[List[float]] = [
        [float(v) for v in df[target].dropna().tolist()] for df in all_dfs
    ]

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
