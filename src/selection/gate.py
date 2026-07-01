"""gate — keep designs that satisfy a set of boolean conditions.

A gate is a plausibility filter: the rows that pass are the point; finer differences are
noise or unreliable (see the funnel's gate-vs-rank principle). Conditions are ANDed at the
top level; each condition is either a LEAF (``{"col": ..., "<op>": value}``) or a GROUP
(``{"all_of": [...]}`` / ``{"any_of": [...]}``) nesting further conditions.

Leaf operators (one per condition): eq, ne, lt, le, gt, ge, in (value = list), not_in,
notnull, isnull, between (value = [lo, hi], inclusive). A MISSING value (NaN) fails every
leaf except isnull — a design lacking a gated metric cannot pass it.

A condition may carry a ``when`` key (itself a condition): rows NOT matching ``when``
auto-pass, rows matching must satisfy the rest — used for class-specific filters (e.g. a
c10-only geometry gate that does not apply to c0/c1).

``apply_gate`` adds a boolean ``gate_pass`` column and (by default) returns only the passing
rows, plus a report dict with the in/out counts and per-condition pass counts for the
provenance manifest.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import pandas as pd

_LEAF_OPS = {"eq", "ne", "lt", "le", "gt", "ge", "in", "not_in", "notnull", "isnull", "between"}


def _describe(cond: dict) -> str:
    if "when" in cond:
        body = {k: v for k, v in cond.items() if k != "when"}
        return f"IF {_describe(cond['when'])} THEN {_describe(body)}"
    if "all_of" in cond:
        return "(" + " AND ".join(_describe(c) for c in cond["all_of"]) + ")"
    if "any_of" in cond:
        return "(" + " OR ".join(_describe(c) for c in cond["any_of"]) + ")"
    col = cond.get("col", "?")
    for op in _LEAF_OPS:
        if op in cond:
            return f"{col} {op} {cond[op]}" if op not in ("notnull", "isnull") else f"{col} {op}"
    return f"{col} <malformed>"


def _eval_leaf(df: pd.DataFrame, cond: dict) -> pd.Series:
    col = cond.get("col")
    if col is None:
        raise ValueError(f"gate leaf missing 'col': {cond}")
    if col not in df.columns:
        raise ValueError(f"gate condition references unknown column '{col}'. "
                         f"available: {list(df.columns)}")
    s = df[col]
    present = s.notna()
    if "notnull" in cond:
        return present
    if "isnull" in cond:
        return ~present
    if "eq" in cond:
        target = cond["eq"]
        # Booleans in a merged CSV may arrive as the strings "True"/"False".
        if isinstance(target, bool):
            norm = s.map(lambda v: str(v).strip().lower() in ("true", "1", "1.0")
                         if pd.notna(v) else False)
            return norm == target
        return present & (s == target)
    if "ne" in cond:
        return present & (s != cond["ne"])
    if "in" in cond:
        return present & s.isin(cond["in"])
    if "not_in" in cond:
        return present & ~s.isin(cond["not_in"])
    if "between" in cond:
        lo, hi = cond["between"]
        return present & (s.astype(float) >= lo) & (s.astype(float) <= hi)
    num = pd.to_numeric(s, errors="coerce")
    if "lt" in cond:
        return num.notna() & (num < cond["lt"])
    if "le" in cond:
        return num.notna() & (num <= cond["le"])
    if "gt" in cond:
        return num.notna() & (num > cond["gt"])
    if "ge" in cond:
        return num.notna() & (num >= cond["ge"])
    raise ValueError(f"gate leaf has no recognised operator ({_LEAF_OPS}): {cond}")


def _eval_condition(df: pd.DataFrame, cond: dict) -> pd.Series:
    if "when" in cond:
        # Conditional: rows NOT matching `when` auto-pass; rows matching must satisfy the
        # rest of the condition. Used for class-specific filters (e.g. c10-only geometry).
        when_mask = _eval_condition(df, cond["when"])
        body = {k: v for k, v in cond.items() if k != "when"}
        return (~when_mask) | _eval_condition(df, body)
    if "all_of" in cond:
        mask = pd.Series(True, index=df.index)
        for c in cond["all_of"]:
            mask &= _eval_condition(df, c)
        return mask
    if "any_of" in cond:
        mask = pd.Series(False, index=df.index)
        for c in cond["any_of"]:
            mask |= _eval_condition(df, c)
        return mask
    return _eval_leaf(df, cond)


def apply_gate(df: pd.DataFrame, conditions: List[dict],
               keep_only_passing: bool = True) -> Tuple[pd.DataFrame, Dict]:
    """Add ``gate_pass`` (top-level AND of ``conditions``) and optionally drop failers."""
    mask = pd.Series(True, index=df.index)
    per_condition = []
    for cond in conditions:
        cmask = _eval_condition(df, cond)
        per_condition.append({"condition": _describe(cond), "passed": int(cmask.sum())})
        mask &= cmask
    out = df.copy()
    out["gate_pass"] = mask.values
    report = {"op": "gate", "n_in": len(df), "n_pass": int(mask.sum()),
              "conditions": per_condition}
    if keep_only_passing:
        out = out[out["gate_pass"]].drop(columns=["gate_pass"])
    return out, report
