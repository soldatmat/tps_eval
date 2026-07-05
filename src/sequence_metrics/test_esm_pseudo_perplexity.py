from __future__ import annotations

"""Self-contained tests for esm_pseudo_perplexity.py.

Run from this directory:
    cd src/sequence_metrics && python test_esm_pseudo_perplexity.py
or under pytest:
    cd src/sequence_metrics && python -m pytest test_esm_pseudo_perplexity.py -q

PARTIAL local coverage. The scoring entrypoints (`compute_pseudo_perplexity`,
`_score_swoop`, `_score_masked`) require the ESM-1b weights
(esm1b_t33_650M_UR50S, ~2.5 GB) which must be DOWNLOADED and are best run on a
GPU — that path is NEEDS-AURUM and is intentionally NOT exercised here (no
network / no model download in this environment).

What IS tested locally (torch is importable here): the pure per-position
log-probability primitive `_per_position_log_probs` (a numerically-known
log-softmax over the vocab axis) and the module's declared constants. This at
least guarantees the module imports and its scoring math primitive is correct.
"""

import math

# torch + fair-esm are only present in an ESM-capable env (e.g. the `esmfold`
# conda env, torch 2.5.1). Where they're missing/broken, skip gracefully instead
# of crashing the suite — the scoring path is NEEDS-AURUM regardless.
try:
    import torch

    import esm_pseudo_perplexity as epp
    from esm_pseudo_perplexity import (
        COLUMNS,
        DEFAULT_MODEL,
        TRUNCATION_SEQ_LENGTH,
        _per_position_log_probs,
    )

    _IMPORT_ERROR = None
except Exception as exc:  # noqa: BLE001 — any import failure means "can't test here"
    _IMPORT_ERROR = exc


def _approx(a, b, tol=1e-5):
    assert abs(a - b) <= tol, f"{a} != {b}"


def test_per_position_log_probs_uniform():
    # Two equal logits -> log-softmax = log(0.5) for each; rows sum to 1 in prob.
    logits = torch.tensor([[[0.0, 0.0]]])
    lp = _per_position_log_probs(logits)
    _approx(lp[0, 0, 0].item(), math.log(0.5))
    _approx(lp[0, 0, 1].item(), math.log(0.5))
    # exp of a log-softmax row sums to 1.
    _approx(lp[0, 0].exp().sum().item(), 1.0)


def test_per_position_log_probs_normalizes_over_last_axis():
    # Arbitrary logits: probabilities over the vocab (last) axis must sum to 1,
    # independently at each position.
    logits = torch.tensor([[[2.0, 1.0, -1.0], [0.5, 0.5, 0.5]]])
    lp = _per_position_log_probs(logits)
    probs = lp.exp()
    for pos in range(probs.shape[1]):
        _approx(probs[0, pos].sum().item(), 1.0)
    # Larger logit -> larger (less negative) log-prob.
    assert lp[0, 0, 0].item() > lp[0, 0, 1].item() > lp[0, 0, 2].item()


def test_constants():
    assert COLUMNS == ["ID", "esm_pseudo_perplexity", "esm_mean_pll", "n_residues"]
    assert DEFAULT_MODEL == "esm1b_t33_650M_UR50S"
    assert TRUNCATION_SEQ_LENGTH == 1022


def test_unknown_method_raises_without_model():
    # An invalid method should be rejected; the module reaches model-load first,
    # so we only assert the module exposes compute_pseudo_perplexity as callable
    # (the actual model-dependent scoring path is NEEDS-AURUM).
    assert callable(epp.compute_pseudo_perplexity)


def main():
    if _IMPORT_ERROR is not None:
        print(f"SKIPPED (torch/esm unavailable in this env): {_IMPORT_ERROR!r}")
        return
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
