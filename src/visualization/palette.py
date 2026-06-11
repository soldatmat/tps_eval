"""First-cyclization-class colour palette, grouped by substrate type.

Each substrate type is one hue family, ordered by precursor carbon number
(mono C10 -> sesqui C15 -> di C20 -> sester C25 -> sterol/triterpene C30);
shades within a family distinguish the classes. Shared by every landscape map
so the ESM / SaProt / EE / similarity maps are directly comparable.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

N_CLASSES = 22

# substrate type per first-cyclization class id (MARTS-DB 2026-04-12, 22-class)
SUBSTRATE: dict[int, str] = {}
for _c in (0, 1, 2, 4, 10, 19):
    SUBSTRATE[_c] = "sesqui"
for _c in (5, 11, 20):
    SUBSTRATE[_c] = "mono"
for _c in (3, 6, 7, 8, 9, 13, 14, 15, 16):
    SUBSTRATE[_c] = "di"
for _c in (17, 18, 21):
    SUBSTRATE[_c] = "sester"
for _c in (12,):
    SUBSTRATE[_c] = "sterol"

# legend / build order = ascending precursor carbon count
TYPE_ORDER = ["mono", "sesqui", "di", "sester", "sterol"]
TYPE_CMAP = {"sesqui": "Blues", "mono": "Greens", "di": "Reds", "sester": "Purples"}
TYPE_FIXED = {"sterol": (0.40, 0.26, 0.13)}  # brown singleton (triterpene/OSC)


def build_class_colors():
    """Return ({class_id: rgba}, {type: [class_ids]}), shaded within each family."""
    by_type = {t: sorted(c for c, tt in SUBSTRATE.items() if tt == t) for t in TYPE_ORDER}
    colors = {}
    for t in TYPE_ORDER:
        ids = by_type[t]
        if t in TYPE_FIXED:
            for cid in ids:
                colors[cid] = (*TYPE_FIXED[t], 1.0)
            continue
        cmap = plt.get_cmap(TYPE_CMAP[t])
        k = len(ids)
        shades = [0.70] if k == 1 else list(np.linspace(0.42, 0.92, k))
        for cid, s in zip(ids, shades):
            colors[cid] = cmap(s)
    return colors, by_type


def make_palette(n: int = N_CLASSES) -> ListedColormap:
    """ListedColormap over class ids 0..n-1, grouped by substrate type."""
    colors, _ = build_class_colors()
    return ListedColormap([colors[c] for c in range(n)])


def grouped_legend(fig, cmap=None, **legend_kw):
    """Attach a class legend grouped by substrate type, parked outside right."""
    if cmap is None:
        cmap = make_palette()
    _, by_type = build_class_colors()
    handles, labels = [], []
    for t in TYPE_ORDER:
        ids = by_type[t]
        handles.append(Line2D([0], [0], linestyle="", marker="", alpha=0))
        labels.append(f"$\\bf{{{t}}}$  ({len(ids)} cls)")
        for c in ids:
            handles.append(Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor=cmap(c), markeredgecolor="none",
                                  markersize=7))
            labels.append(f"   {c:>2d}")
    kw = dict(loc="center left", bbox_to_anchor=(1.0, 0.5),
              bbox_transform=fig.transFigure, fontsize=8, title_fontsize=9,
              ncol=1, frameon=True, framealpha=0.9, handletextpad=0.4,
              labelspacing=0.35)
    kw.update(legend_kw)
    fig.legend(handles=handles, labels=labels,
               title="first-cyclization class\n(grouped by substrate type)", **kw)
