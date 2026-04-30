from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.load_results import load_results  # noqa: E402
from plot.constants import (  # noqa: E402
    LOAD,
    MAX_VAL,
    MIN_VAL,
    OFFSET,
    THRESHOLD,
    TICKS,
)


def _kde_curve(values: List[float], grid: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.zeros_like(grid)
    if arr.size == 1 or np.allclose(arr, arr[0]):
        # gaussian_kde fails on zero variance; emit a narrow spike at the value.
        spike = np.zeros_like(grid)
        idx = int(np.argmin(np.abs(grid - arr[0])))
        spike[idx] = 1.0
        return spike
    kde = gaussian_kde(arr)
    return kde(grid)


def density_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    target: str,
    *,
    direction: str = "vertical",
    save_dir: Optional[str] = None,
) -> None:
    """Ridge-style density plot of `target` across datasets.

    `direction='vertical'` lays each dataset's KDE along the y-axis, offset on
    x. `direction='horizontal'` does the inverse.
    """
    if direction not in ("vertical", "horizontal"):
        raise ValueError(f"direction must be 'vertical' or 'horizontal', got {direction!r}")

    load_list = LOAD[target]
    min_val = MIN_VAL[target]
    max_val = MAX_VAL[target]
    ticks = TICKS[target]
    threshold = THRESHOLD[target]
    offset = OFFSET[target]

    all_dfs = [load_results(fp, load=load_list) for fp in fasta_paths]
    all_data: List[List[float]] = [
        [float(v) for v in df[target].dropna().tolist()] for df in all_dfs
    ]

    grid = np.linspace(min_val, max_val, 512)

    fig, ax = plt.subplots()
    n = len(all_data)
    label_positions = [(i + 1) * offset for i in range(n)]

    for i in range(n - 1, -1, -1):
        density = _kde_curve(all_data[i], grid)
        baseline = (i + 1) * offset
        if direction == "vertical":
            ax.fill_betweenx(
                grid,
                baseline,
                baseline + density,
                color=data_colors[i],
                edgecolor="black",
                linewidth=1,
            )
        else:
            ax.fill_between(
                grid,
                baseline,
                baseline + density,
                color=data_colors[i],
                edgecolor="black",
                linewidth=1,
            )

    ax.set_title(target + " density")

    if direction == "vertical":
        ax.set_yticks(ticks)
        ax.set_xticks(label_positions)
        ax.set_xticklabels(data_names)
        ax.set_ylim(min_val, max_val)
        ax.yaxis.grid(True, linestyle="--", color="gray")
        ax.set_axisbelow(True)
        if threshold is not None:
            ax.axhline(threshold, color="red", linestyle="--", linewidth=2)
    else:
        ax.set_xticks(ticks)
        ax.set_yticks(label_positions)
        ax.set_yticklabels(data_names)
        ax.set_xlim(min_val, max_val)
        ax.xaxis.grid(True, linestyle="--", color="gray")
        ax.set_axisbelow(True)
        if threshold is not None:
            ax.axvline(threshold, color="red", linestyle="--", linewidth=2)

    fig.tight_layout()

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, target + "_density_" + direction + ".png"))

    plt.close(fig)
