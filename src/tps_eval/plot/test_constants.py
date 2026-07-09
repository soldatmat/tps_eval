from __future__ import annotations

"""Self-contained tests for plot/constants.py (the plot-target metadata tables).

Run from src/ so the package-style imports resolve like the plot runner does:
    cd src && python -m plot.test_constants
or under pytest:
    cd src && python -m pytest plot/test_constants.py -q

Pure metadata-consistency checks (no plotting): the numeric-target tables agree
with each other (every keyed target is a known TARGET; MIN<MAX; ticks within
bounds; thresholds within bounds), and the structure-branch column tables are
internally consistent (a suffix's numeric and categorical column sets are
disjoint, no dupes). Also locks in the domain-composition column-name fix
(the CSV emits n_<type> counts, not bare type names).
"""

import sys
from pathlib import Path


import numpy as np

from tps_eval.plot import constants as C


def test_targets_unique():
    assert len(C.TARGETS) == len(set(C.TARGETS)), "duplicate entries in TARGETS"


def test_every_target_is_loadable():
    """Each numeric TARGET must map to a LOAD tool CSV (else it can never load)."""
    missing = [t for t in C.TARGETS if t not in C.LOAD]
    assert not missing, f"TARGETS with no LOAD entry: {missing}"


def test_min_less_than_max():
    for t in C.MIN_VAL:
        assert t in C.MAX_VAL, f"{t} in MIN_VAL but not MAX_VAL"
        assert C.MIN_VAL[t] < C.MAX_VAL[t], f"{t}: min >= max"
    for t in C.MAX_VAL:
        assert t in C.MIN_VAL, f"{t} in MAX_VAL but not MIN_VAL"


def test_ticks_within_bounds():
    """Every explicit tick array stays inside that target's [min, max] window."""
    for t, ticks in C.TICKS.items():
        arr = np.asarray(ticks, dtype=float)
        assert arr.ndim == 1 and arr.size > 0, t
        # monotonically increasing
        assert np.all(np.diff(arr) > 0), f"{t}: ticks not strictly increasing"
        if t in C.MIN_VAL:
            assert arr.min() >= C.MIN_VAL[t] - 1e-9, f"{t}: tick below min"
            assert arr.max() <= C.MAX_VAL[t] + 1e-9, f"{t}: tick above max"


def test_threshold_within_bounds():
    for t, thr in C.THRESHOLD.items():
        if thr is None:
            continue
        if t in C.MIN_VAL:
            assert C.MIN_VAL[t] <= thr <= C.MAX_VAL[t], f"{t}: threshold outside axis"


def test_structure_numeric_categorical_disjoint_per_suffix():
    """A column must not be both a numeric and a categorical target for a suffix."""
    for suffix, num_cols in C.STRUCTURE_NUMERIC.items():
        assert len(num_cols) == len(set(num_cols)), f"{suffix}: dup numeric cols"
        cat_cols = C.STRUCTURE_CATEGORICAL.get(suffix, [])
        assert len(cat_cols) == len(set(cat_cols)), f"{suffix}: dup categorical cols"
        overlap = set(num_cols) & set(cat_cols)
        assert not overlap, f"{suffix}: columns both numeric and categorical: {overlap}"


def test_structure_suffixes_look_like_csv_suffixes():
    for suffix in list(C.STRUCTURE_NUMERIC) + list(C.STRUCTURE_CATEGORICAL):
        assert suffix.startswith("_") and suffix.endswith(".csv"), suffix


def test_domain_composition_uses_prefixed_count_columns():
    """Regression: the domain_composition CSV emits n_<type> counts (n_alpha ...
    n_zeta), not bare type names. A mismatch silently skips every count plot."""
    cols = C.STRUCTURE_NUMERIC["_domain_composition.csv"]
    expected = {
        "n_domains", "n_alpha", "n_beta", "n_gamma",
        "n_ids", "n_delta", "n_epsilon", "n_zeta",
    }
    assert set(cols) == expected, cols
    # None of the stale bare names survive.
    assert not ({"alpha", "beta", "gamma", "delta", "epsilon", "ids",
                 "terpene_synth_C"} & set(cols))


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
