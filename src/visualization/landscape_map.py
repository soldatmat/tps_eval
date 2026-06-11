"""Render first-cyclization-class-coloured 2D landscape maps.

Given precomputed 2D coordinates (from dimensionality_reduction) and a per-row
class label, draw one scatter panel per layout, coloured with the shared
substrate-type palette and a grouped legend parked outside the data area.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .palette import make_palette, grouped_legend, N_CLASSES


def render_panels(panels, y, suptitle, out_png, footnote=None, figsize=None,
                  s=10, alpha=0.8):
    """Render a multi-panel class-coloured landscape figure.

    panels   : list of (coords (n,2), panel_title)
    y        : (n,) int first-cyclization class per row
    out_png  : output path (PNG)
    """
    cmap = make_palette()
    k = len(panels)
    fig, axes = plt.subplots(1, k, figsize=figsize or (6.6 * k, 6.0),
                             constrained_layout=True)
    if k == 1:
        axes = [axes]
    for ax, (coords, title) in zip(axes, panels):
        ax.scatter(coords[:, 0], coords[:, 1], c=y, cmap=cmap,
                   vmin=-0.5, vmax=N_CLASSES - 0.5, s=s, alpha=alpha, linewidths=0)
        ax.set_title(title, fontsize=11)
        ax.tick_params(labelsize=8)
    grouped_legend(fig, cmap)
    fig.suptitle(suptitle, fontsize=12.5, fontweight="bold")
    if footnote:
        fig.text(0.01, 0.005, footnote, fontsize=9, style="italic",
                 ha="left", va="bottom")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_png
