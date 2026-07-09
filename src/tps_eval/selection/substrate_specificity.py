"""substrate_specificity — keep designs that are ON-target AND specific for one substrate.

A campaign targets one substrate (the prenyl-diphosphate the enzyme is being designed to
act on). This selection op reads the EnzymeExplorer sequence-only per-substrate scores and
keeps a design iff:

    EE[target] >= t_hi                 (on-target: confidently the target substrate)
    AND  EE[other] <= ceiling  for every other scored substrate   (off-target: not also
                                       strongly predicted for a competing substrate)

The off-target ceiling is deliberately RELAXED (lenient/high) — it prunes designs that are
nearly as good on a competing substrate, not designs with any trace off-target signal.

EE score columns come in two formats and both are handled:
  * structured output — "<SMILES> (<name>)" columns (e.g. "... (Farnesyl pyrophosphate)"),
    matched by the exact parenthetical name (see _EE_NAME_TO_CODE);
  * console `predict_sequences_only` output — "<CODE>_score" columns (FPP_score, GPP_score, ...).
EE-specific codes are folded onto the shared substrate vocabulary (CPP -> GGPP, 2xFPP -> EDSQ),
matching src/knn/substrate_class.py, so neither the copalyl nor 2xFPP column is lost. When a
substrate maps to several columns (e.g. GGPP + folded CPP) the per-substrate score is their max.

Defaults (t_hi=0.5, off-target ceiling=0.35): EE seq-only scores are softmax-like
probabilities over ~9-11 substrate classes summing to ~1, so an on-target design puts the
majority of the mass on the target (>=0.5) while every competing substrate stays well below
(a typical runner-up is ~0.15-0.20, so a 0.35 ceiling is lenient). Both are tunable per spec.

Missing data:
  * If EE does not score the target substrate AT ALL (no column for it — e.g. DMAPP/C35/IDS,
    which EE seq-only does not emit) the gate cannot be evaluated -> ValueError (spec/config
    error, fail loud). Use a plain `gate` instead for those substrates.
  * A per-ROW NaN target score -> that design fails the on-target test (dropped), per design.
  * A per-ROW NaN off-target score -> treated as not-exceeding (passes that ceiling).
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

DEFAULT_T_HI = 0.5
DEFAULT_T_OFF = 0.35

_EE_SCORE_SUFFIX = "_score"
_EE_NON_SUBSTRATE = {"TPS", "isTPS"}
# Fold EE-specific codes onto the shared substrate-label vocabulary (mirrors substrate_class).
_EE_CLASS_FOLD = {
    "CPP": "GGPP",     # copalyl-PP -> C20 diterpene
    "2xFPP": "EDSQ",   # 2xFPP -> squalene / 2,3-epoxysqualene (C30)
}
# Exact parenthetical name (lowercased) -> substrate code, for the structured EE output.
_EE_NAME_TO_CODE = {
    "dimethylallyl pyrophosphate": "DMAPP",
    "geranyl pyrophosphate": "GPP",
    "farnesyl pyrophosphate": "FPP",
    "geranylgeranyl pyrophosphate": "GGPP",
    "geranylfarnesyl pyrophosphate": "GFPP",
    "(s)-2,3-epoxysqualene": "EDSQ",
    "2,3-epoxysqualene": "EDSQ",
    "copalyl diphosphate": "GGPP",
    "2x farnesyl pyrophosphate": "EDSQ",
    "2x geranylgeranyl pyrophosphate": "2xGGPP",
}

_PAREN_NAME = re.compile(r"\(([^)]*)\)\s*$")


def _fold(code: str) -> str:
    return _EE_CLASS_FOLD.get(code, code)


def ee_columns_by_substrate(df: pd.DataFrame) -> Dict[str, List[str]]:
    """Map substrate code -> list of ORIGINAL EE score column names present in ``df``.

    Recognises both the "<CODE>_score" console columns and the "<SMILES> (<name>)"
    structured columns, folding CPP->GGPP and 2xFPP->EDSQ. Column names are matched with
    surrounding whitespace stripped (EE output has trailing spaces on some headers), but the
    ORIGINAL name is returned so the caller can index ``df`` directly.
    """
    out: Dict[str, List[str]] = {}
    for original in df.columns:
        stripped = str(original).strip()
        code: Optional[str] = None
        if stripped.endswith(_EE_SCORE_SUFFIX):
            raw = stripped[: -len(_EE_SCORE_SUFFIX)].strip()
            if raw and raw not in _EE_NON_SUBSTRATE:
                code = raw
        else:
            m = _PAREN_NAME.search(stripped)
            if m:
                code = _EE_NAME_TO_CODE.get(m.group(1).strip().lower())
        if code is None:
            continue
        out.setdefault(_fold(code), []).append(original)
    return out


def _max_over(df: pd.DataFrame, columns: List[str]) -> pd.Series:
    """Row-wise max over the given columns (numeric-coerced), NaN where all are NaN."""
    numeric = pd.concat([pd.to_numeric(df[c], errors="coerce") for c in columns], axis=1)
    return numeric.max(axis=1, skipna=True)


def apply_substrate_specificity(
    df: pd.DataFrame,
    target_substrate: str,
    *,
    t_hi: float = DEFAULT_T_HI,
    t_off: float = DEFAULT_T_OFF,
    per_substrate_ceilings: Optional[Dict[str, float]] = None,
    keep_only_passing: bool = True,
) -> Tuple[pd.DataFrame, Dict]:
    """Keep designs on-target for ``target_substrate`` and specific against the rest.

    Adds a boolean ``specificity_pass`` column and (by default) drops the failers. Returns
    (df, report). Raises ValueError if EE does not score the target substrate at all.
    """
    target = _fold(str(target_substrate).strip().upper())
    ceilings = {str(k).strip().upper(): float(v) for k, v in (per_substrate_ceilings or {}).items()}

    cols_by_code = ee_columns_by_substrate(df)
    if target not in cols_by_code:
        raise ValueError(
            f"substrate_specificity: EE does not score target substrate {target!r}. "
            f"Available EE substrate columns: {sorted(cols_by_code)}. "
            f"(EE seq-only does not emit DMAPP/C35/IDS — use a plain gate for those.)"
        )

    on = _max_over(df, cols_by_code[target])
    on_pass = on.notna() & (on >= t_hi)

    off_pass = pd.Series(True, index=df.index)
    per_offtarget: List[dict] = []
    for code, columns in sorted(cols_by_code.items()):
        if code == target:
            continue
        ceiling = ceilings.get(code, t_off)
        off = _max_over(df, columns)
        # A missing off-target score does not exceed the ceiling -> passes.
        this_pass = off.isna() | (off <= ceiling)
        off_pass &= this_pass
        per_offtarget.append({"substrate": code, "ceiling": ceiling,
                              "passed": int(this_pass.sum())})

    mask = on_pass & off_pass
    out = df.copy()
    out["specificity_pass"] = mask.values

    report = {
        "op": "substrate_specificity",
        "target": target,
        "t_hi": t_hi,
        "t_off": t_off,
        "n_in": len(df),
        "n_pass": int(mask.sum()),
        "n_on_target_pass": int(on_pass.sum()),
        "n_missing_target_score": int(on.isna().sum()),
        "off_target": per_offtarget,
    }
    if keep_only_passing:
        out = out[out["specificity_pass"]].drop(columns=["specificity_pass"])
    return out, report
