"""Self-contained tests for the substrate_specificity selection op.

Run from this directory (flat-module imports resolve like the runners do):
    cd src/selection && python test_substrate_specificity.py
or under pytest:
    cd src/selection && python -m pytest test_substrate_specificity.py -q
"""
import os
import sys

import numpy as np
import pandas as pd


from tps_eval.selection.substrate_specificity import (
    apply_substrate_specificity,
    ee_columns_by_substrate,
)


def _ee_df_structured():
    """EE 'structured' output: '<SMILES> (<name>)' columns (the committed sample format)."""
    return pd.DataFrame({
        "ID": ["a", "b", "c", "d"],
        "smi1 (Farnesyl pyrophosphate)":       [0.60, 0.40, 0.55, 0.52],  # FPP (target)
        "smi2 (Geranyl pyrophosphate)":        [0.10, 0.10, 0.10, 0.10],  # GPP
        "smi3 (Geranylgeranyl pyrophosphate)": [0.10, 0.10, 0.40, 0.10],  # GGPP
        "smi4 (copalyl diphosphate)":          [0.05, 0.05, 0.05, 0.50],  # CPP -> folds to GGPP
        "isTPS ": [0.7, 0.7, 0.7, 0.7],
    })


def test_column_resolution_both_formats_and_folding():
    cols = ee_columns_by_substrate(_ee_df_structured())
    assert set(cols) == {"FPP", "GPP", "GGPP"}, cols
    # copalyl folds into GGPP -> GGPP collects two columns.
    assert len(cols["GGPP"]) == 2, cols["GGPP"]

    score_fmt = pd.DataFrame({
        "ID": ["a"],
        "FPP_score": [0.6], "GPP_score": [0.2],
        "CPP_score": [0.5],          # -> GGPP
        "2xFPP_score": [0.1],        # -> EDSQ
        "TPS_score": [0.9], "isTPS": [0.9],
    })
    cols2 = ee_columns_by_substrate(score_fmt)
    assert set(cols2) == {"FPP", "GPP", "GGPP", "EDSQ"}, cols2
    assert "TPS" not in cols2 and "isTPS" not in cols2
    print("ok column resolution (both formats + CPP/2xFPP folding, TPS/isTPS excluded)")


def test_on_and_off_target_gating_with_folding():
    out, rep = apply_substrate_specificity(_ee_df_structured(), "FPP", t_hi=0.5, t_off=0.35)
    # a: on 0.60, off max 0.10 -> keep. b: on 0.40 < 0.5 -> drop.
    # c: on 0.55 but GGPP 0.40 > 0.35 -> drop. d: on 0.52 but folded CPP 0.50 -> GGPP -> drop.
    assert set(out["ID"]) == {"a"}, set(out["ID"])
    assert rep["n_pass"] == 1 and rep["n_on_target_pass"] == 3
    print("ok on/off-target gating (with CPP fold pushing d over the ceiling)")


def test_missing_target_score_drops_row():
    df = _ee_df_structured()
    df.loc[df["ID"] == "a", "smi1 (Farnesyl pyrophosphate)"] = np.nan
    out, rep = apply_substrate_specificity(df, "FPP", t_hi=0.5, t_off=0.35)
    assert "a" not in set(out["ID"]), set(out["ID"])
    assert rep["n_missing_target_score"] == 1
    print("ok missing target score -> row dropped")


def test_nan_offtarget_passes_ceiling():
    df = pd.DataFrame({
        "ID": ["a"],
        "smi1 (Farnesyl pyrophosphate)": [0.9],
        "smi2 (Geranyl pyrophosphate)": [np.nan],  # NaN off-target -> not exceeding
    })
    out, _ = apply_substrate_specificity(df, "FPP", t_hi=0.5, t_off=0.35)
    assert set(out["ID"]) == {"a"}
    print("ok NaN off-target treated as not-exceeding")


def test_per_substrate_ceiling_override():
    df = _ee_df_structured()
    # Relax the GGPP ceiling to 0.6 -> c (GGPP 0.40) and d (folded 0.50) now pass off-target;
    # both are on-target (0.55, 0.52) -> kept alongside a.
    out, _ = apply_substrate_specificity(
        df, "FPP", t_hi=0.5, t_off=0.35, per_substrate_ceilings={"GGPP": 0.6})
    assert set(out["ID"]) == {"a", "c", "d"}, set(out["ID"])
    print("ok per-substrate ceiling override")


def test_target_not_scored_raises():
    try:
        apply_substrate_specificity(_ee_df_structured(), "DMAPP")
    except ValueError as e:
        assert "DMAPP" in str(e)
        print("ok target not scored by EE -> ValueError")
        return
    raise AssertionError("expected ValueError for unscored target substrate")


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all substrate_specificity tests passed")


if __name__ == "__main__":
    main()
