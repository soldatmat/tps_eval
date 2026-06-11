from __future__ import annotations

import glob
import os
import sys
import traceback
from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib.colors as mcolors

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.load_results import load_results, load_structure_results  # noqa: E402
from plot.boxplot_comparison import boxplot_comparison  # noqa: E402
from plot.categorical_comparison import categorical_comparison  # noqa: E402
from plot.constants import (  # noqa: E402
    CATEGORICAL_TARGETS,
    STRUCTURE_CATEGORICAL,
    STRUCTURE_NUMERIC,
    TARGETS,
)
from plot.density_comparison import density_comparison  # noqa: E402

# Fallback palette for color names matplotlib can't resolve at all.
_FALLBACK_COLORS = ("dodgerblue", "goldenrod", "seagreen", "tomato", "purple", "gray")

# Columns in *_motifs.csv that are bookkeeping, not motif-presence booleans.
_MOTIF_NON_TARGET_COLS = {"ID", "sequence", "msa"}


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


def _plot_numeric(
    fasta_paths, data_names, data_colors, target, *, dfs=None, save_dir=None
) -> None:
    """Boxplot + vertical-density pair for one numeric target."""
    print(f"Generating boxplot for target: {target}...")
    boxplot_comparison(
        fasta_paths, data_names, data_colors, target, dfs=dfs, save_dir=save_dir
    )
    print(f"Generating density plot for target: {target}...")
    density_comparison(
        fasta_paths,
        data_names,
        data_colors,
        target,
        direction="vertical",
        dfs=dfs,
        save_dir=save_dir,
    )


def _discover_motif_targets(fasta_paths: Sequence[str]) -> List[str]:
    """Boolean motif-presence columns from the *_motifs.csv files.

    The motif CSV's column names ARE the regex patterns, so the set of targets
    is data-dependent. We union the boolean columns across datasets. Returns []
    (and is a no-op for plotting) when no motifs CSV is present.
    """
    cols: List[str] = []
    for fp in fasta_paths:
        try:
            df = load_results(fp, load=["motifs"])
        except FileNotFoundError:
            continue
        for c in df.columns:
            if c in _MOTIF_NON_TARGET_COLS or c in cols:
                continue
            # Motif columns are booleans; guard against any stray numeric columns.
            if df[c].dropna().isin([True, False, 0, 1]).all():
                cols.append(c)
    return cols


def plot_comparison(
    fasta_paths: Sequence[str],
    data_names: Sequence[str],
    data_colors: Sequence[str],
    *,
    targets: Optional[Sequence[str]] = None,
    save_dir: Optional[str] = None,
) -> None:
    data_colors = [
        _normalize_color(c, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
        for i, c in enumerate(data_colors)
    ]

    # ------------------------------------------------------------------ #
    # 1. Sequence-branch NUMERIC targets (comparison across all datasets)
    # ------------------------------------------------------------------ #
    numeric_targets = TARGETS if targets is None else list(targets)
    for target in numeric_targets:
        try:
            _plot_numeric(
                fasta_paths, data_names, data_colors, target, save_dir=save_dir
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

    # When the caller hand-picked `--targets`, honor that exactly and skip the
    # auto-discovered categorical / structure passes below.
    if targets is not None:
        return

    # ------------------------------------------------------------------ #
    # 2. Sequence-branch CATEGORICAL / BOOLEAN targets (count plots)
    # ------------------------------------------------------------------ #
    categorical_targets = list(CATEGORICAL_TARGETS.items())
    # Dynamic motif-presence boolean columns (all live in the *_motifs.csv).
    for motif_col in _discover_motif_targets(fasta_paths):
        categorical_targets.append((motif_col, ["motifs"]))

    for target, load_list in categorical_targets:
        try:
            print(f"Generating count plot for categorical target: {target}...")
            categorical_comparison(
                fasta_paths,
                data_names,
                data_colors,
                target,
                load_list=load_list,
                save_dir=save_dir,
            )
        except FileNotFoundError as e:
            print(f"  [skip] target {target}: missing input ({e.filename or e})")
        except Exception as e:
            print(f"Error while plotting for categorical target {target}: {e}")
            print("Stacktrace:")
            traceback.print_exc()

    # ------------------------------------------------------------------ #
    # 3. STRUCTURE-branch targets (generated set only → single distribution)
    # ------------------------------------------------------------------ #
    # Structures exist for the generated set only; their CSVs are named
    # `<structs_dir>_<suffix>.csv` and live in the input directory (a sibling of
    # the fasta files). We attach them to the generated dataset, which the
    # orchestrator passes LAST (data_names = "train", "generated").
    gen_fasta = fasta_paths[-1]
    gen_name = data_names[-1]
    gen_color = data_colors[-1]
    input_dir = os.path.dirname(os.path.abspath(gen_fasta))
    _plot_structure_targets(input_dir, gen_name, gen_color, save_dir=save_dir)


def _find_structure_csv(input_dir: str, suffix: str) -> Optional[str]:
    """Locate the single `*<suffix>` structure CSV in `input_dir`, or None."""
    matches = sorted(glob.glob(os.path.join(input_dir, "*" + suffix)))
    return matches[0] if matches else None


def _plot_structure_targets(
    input_dir: str, gen_name: str, gen_color: str, *, save_dir: Optional[str]
) -> None:
    # NUMERIC structure metrics → single-distribution boxplot + density.
    for suffix, columns in STRUCTURE_NUMERIC.items():
        csv_path = _find_structure_csv(input_dir, suffix)
        if csv_path is None:
            for target in columns:
                print(f"  [skip] target {target}: missing input (*{suffix})")
            continue
        struct_df = load_structure_results(csv_path)
        for target in columns:
            if target not in struct_df.columns:
                print(f"  [skip] target {target}: column absent in {os.path.basename(csv_path)}")
                continue
            try:
                _plot_numeric(
                    [gen_fasta_placeholder()],
                    [gen_name],
                    [gen_color],
                    target,
                    dfs=[struct_df],
                    save_dir=save_dir,
                )
            except Exception as e:
                print(f"Error while plotting for structure target {target}: {e}")
                print("Stacktrace:")
                traceback.print_exc()

    # CATEGORICAL structure metrics → single-dataset count plot.
    for suffix, columns in STRUCTURE_CATEGORICAL.items():
        csv_path = _find_structure_csv(input_dir, suffix)
        if csv_path is None:
            for target in columns:
                print(f"  [skip] target {target}: missing input (*{suffix})")
            continue
        struct_df = load_structure_results(csv_path)
        for target in columns:
            if target not in struct_df.columns:
                print(f"  [skip] target {target}: column absent in {os.path.basename(csv_path)}")
                continue
            try:
                print(f"Generating count plot for categorical target: {target}...")
                categorical_comparison(
                    [gen_fasta_placeholder()],
                    [gen_name],
                    [gen_color],
                    target,
                    dfs=[struct_df],
                    save_dir=save_dir,
                )
            except Exception as e:
                print(f"Error while plotting for categorical target {target}: {e}")
                print("Stacktrace:")
                traceback.print_exc()


def gen_fasta_placeholder() -> str:
    """Placeholder fasta path for structure plots.

    Structure plots pass `dfs=` so the fasta path is never read; this keeps the
    positional signature of the plot helpers happy.
    """
    return "__structure__.fasta"
