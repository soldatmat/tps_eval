from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.colors as mcolors

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from plot.boxplot_comparison import boxplot_comparison  # noqa: E402
from plot.constants import TARGETS  # noqa: E402
from plot.density_comparison import density_comparison  # noqa: E402

# Fallback palette for color names matplotlib can't resolve at all.
_FALLBACK_COLORS = ("dodgerblue", "goldenrod", "seagreen", "tomato", "purple", "gray")


def _normalize_color(color: str, fallback: str) -> str:
    """Map a possibly-invalid color name to a matplotlib-valid one.

    Accepts X11/R-style names with a trailing shade digit (e.g. 'dodgerblue3',
    'goldenrod1') that matplotlib doesn't know by stripping the digit to the base
    name ('dodgerblue', 'goldenrod'). If that still isn't a valid color, returns
    `fallback`.
    """
    if mcolors.is_color_like(color):
        return color
    stripped = color.rstrip("0123456789")
    if stripped and mcolors.is_color_like(stripped):
        print(f"  [plot] normalized color '{color}' -> '{stripped}'")
        return stripped
    print(f"  [plot] color '{color}' is not a valid matplotlib color; using '{fallback}'")
    return fallback


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

    data_colors = [
        _normalize_color(c, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
        for i, c in enumerate(data_colors)
    ]

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
        except FileNotFoundError as e:
            # A required per-target input CSV (e.g. *_soluprot.csv,
            # *_enzyme_explorer.csv) wasn't produced for this run — skip the
            # target cleanly instead of failing the whole plot pass.
            print(f"  [skip] target {target}: missing input ({e.filename or e})")
        except Exception as e:  # mirror Julia's per-target try/catch
            print(f"Error while plotting for target {target}: {e}")
            print("Stacktrace:")
            traceback.print_exc()
