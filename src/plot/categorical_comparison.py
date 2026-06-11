from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.load_results import load_results  # noqa: E402


def _category_label(value) -> str:
    """Render a category value as a stable bar label."""
    if isinstance(value, bool):
        return str(value)
    # pandas may give numpy bools / NA; normalize to a readable string.
    s = str(value).strip()
    return s if s != "" else "(none)"


def categorical_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    target: str,
    *,
    load_list: Optional[Sequence[str]] = None,
    dfs: Optional[Sequence] = None,
    save_dir: Optional[str] = None,
) -> None:
    """Grouped count plot of a categorical/boolean `target` across datasets.

    Each dataset becomes a group of bars; bar height = number of rows with that
    category value. NaN/NA values are dropped. Handles one or many datasets
    (structure-only categoricals plot the single generated series). Pass `dfs`
    to bypass fasta-based loading (structure targets); otherwise `load_list`
    selects which sequence CSVs to merge.
    """
    if dfs is None:
        dfs = [load_results(fp, load=load_list) for fp in fasta_paths]
    all_dfs = dfs

    # Per-dataset value->count, in first-seen order across all datasets.
    per_dataset_counts: List[dict] = []
    categories: List[str] = []
    for df in all_dfs:
        counts: dict = {}
        if target in df.columns:
            series = df[target].dropna()
            for value in series.tolist():
                label = _category_label(value)
                counts[label] = counts.get(label, 0) + 1
                if label not in categories:
                    categories.append(label)
        per_dataset_counts.append(counts)

    if not categories:
        # No data for this target in any dataset → nothing to draw.
        print(f"  [skip] categorical target {target}: no values found")
        return

    categories = sorted(categories)
    n_datasets = len(all_dfs)
    n_cats = len(categories)

    fig, ax = plt.subplots()
    x = np.arange(n_cats)
    total_width = 0.8
    bar_width = total_width / max(n_datasets, 1)

    for i, (name, counts) in enumerate(zip(data_names, per_dataset_counts)):
        heights = [counts.get(cat, 0) for cat in categories]
        offsets = x - total_width / 2 + bar_width * (i + 0.5)
        ax.bar(
            offsets,
            heights,
            width=bar_width,
            color=data_colors[i],
            edgecolor="black",
            label=name,
        )

    ax.set_title(target + " counts")
    ax.set_ylabel("count")
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=45, ha="right")
    ax.yaxis.grid(True, linestyle="--", color="gray")
    ax.set_axisbelow(True)
    if n_datasets > 1:
        ax.legend()

    fig.tight_layout()

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, target + "_counts.png"))

    plt.close(fig)
