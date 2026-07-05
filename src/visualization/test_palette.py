from __future__ import annotations

"""Self-contained tests for visualization/palette.py.

Run from this directory (so the package-relative imports resolve like the
package does):
    cd src && python -m visualization.test_palette
or under pytest:
    cd src && python -m pytest visualization/test_palette.py -q

These lock in the substrate-type colour map's structural invariants: the
class-id -> substrate partition is a complete, disjoint cover of 0..N_CLASSES-1;
every class gets a valid, opaque RGBA colour; and make_palette yields exactly
N_CLASSES distinct colours indexable by class id. Pure data / colour logic only
(the legend + panel renderers are drawing-only and are not exercised).
"""

import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import matplotlib

matplotlib.use("Agg")

from matplotlib.colors import ListedColormap, to_hex

from visualization.palette import (
    N_CLASSES,
    SUBSTRATE,
    TYPE_ORDER,
    build_class_colors,
    make_palette,
)


def test_substrate_partition_is_complete_and_disjoint():
    """Every class id 0..N_CLASSES-1 is assigned to exactly one substrate type."""
    assert set(SUBSTRATE.keys()) == set(range(N_CLASSES))
    assert len(SUBSTRATE) == N_CLASSES
    assert set(SUBSTRATE.values()) <= set(TYPE_ORDER)


def test_build_class_colors_covers_all_ids_disjointly():
    colors, by_type = build_class_colors()
    # by_type partitions the ids: union == 0..N-1, pairwise disjoint.
    assert set(by_type.keys()) == set(TYPE_ORDER)
    all_ids = [cid for ids in by_type.values() for cid in ids]
    assert sorted(all_ids) == list(range(N_CLASSES))  # complete + no duplicates
    assert set(colors.keys()) == set(range(N_CLASSES))
    # Each within-family id list is sorted ascending (build order == class-id order).
    for ids in by_type.values():
        assert ids == sorted(ids)


def test_colors_are_valid_opaque_rgba():
    colors, _ = build_class_colors()
    for cid, rgba in colors.items():
        assert len(rgba) == 4, (cid, rgba)
        for ch in rgba:
            assert 0.0 <= ch <= 1.0, (cid, rgba)
        assert rgba[3] == 1.0, f"class {cid} not opaque"  # alpha == 1


def test_make_palette_len_and_uniqueness():
    cmap = make_palette()
    assert isinstance(cmap, ListedColormap)
    assert cmap.N == N_CLASSES
    # Colour per class id is unique (no two classes collide to the same hex).
    hexes = [to_hex(cmap(c)) for c in range(N_CLASSES)]
    assert len(set(hexes)) == N_CLASSES, "duplicate class colours"


def test_make_palette_subset_n():
    """make_palette(n) with n < N_CLASSES yields exactly n colours."""
    cmap = make_palette(5)
    assert cmap.N == 5


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
