"""Small pandas shims that keep behaviour stable across major versions.

The cluster envs are not upgraded in lockstep (Aurum was on pandas 2.3.3 while a
dev laptop was already on 3.0.3), so a behaviour change that only bites on the
newer major silently becomes "works on the cluster, crashes for whoever upgrades
first". Pin the behaviour here rather than at the call sites.
"""

from __future__ import annotations

import pandas as pd


def group_idxmax_skipna(grouped, column: str) -> pd.Series:
    """``grouped[column].idxmax()``, but an ALL-NA group yields NaN instead of raising.

    pandas < 3 returned NaN for such a group; pandas 3 raises
    ``ValueError: idxmax with skipna=True encountered all NA values in a group``.
    The foldseek best-hit reducers rely on the NaN (a query whose whole score column
    is non-numeric must produce a NaN row, not kill the run), so restore it.

    Implemented by excluding the all-NA groups before the native, vectorized idxmax
    and reindexing them back in as NaN — no per-group Python callback, so this stays
    cheap on the ~300k-group production reductions.

    `grouped` is a DataFrameGroupBy; the result is indexed by group key, matching
    what ``grouped[column].idxmax()`` returns.
    """
    counts = grouped[column].count()  # non-NA count per group
    scorable = counts.index[counts > 0]
    frame = grouped.obj
    keys = grouped.keys
    mask = frame[keys].isin(scorable)
    if not mask.any():
        return pd.Series(pd.NA, index=counts.index, dtype="object")
    idx = frame[mask].groupby(keys)[column].idxmax()
    return idx.reindex(counts.index)
