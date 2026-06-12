"""Tests for extract_pae.confidences_to_npz, focused on the protein-chain
restriction for multi-chain AF3 holo co-folds.

Pure numpy + json + tempfile (no AF3 / GPU). Plain-assert runner matching the
repo convention (e.g. src/structure_metrics/test_interdomain_pae.py): run with
``python test_extract_pae.py`` from this directory.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_pae import confidences_to_npz  # noqa: E402


def _write_job(tmp: str, *, pae, token_res_ids, token_chain_ids, ptm=0.8, iptm=None):
    """Write a fake <job>_confidences.json (+ summary) and return its path."""
    conf = os.path.join(tmp, "job_confidences.json")
    with open(conf, "w") as fh:
        json.dump(
            {
                "pae": np.asarray(pae).tolist(),
                "token_res_ids": list(token_res_ids),
                "token_chain_ids": list(token_chain_ids),
            },
            fh,
        )
    with open(os.path.join(tmp, "job_summary_confidences.json"), "w") as fh:
        json.dump({"ptm": ptm, "iptm": iptm}, fh)
    return conf


def test_holo_multichain_restricted_to_protein_chain():
    """A protein chain (A, 3 tokens) + 3 Mg (B,C,D) + POP (E) -> PAE sliced to the
    3 protein tokens, residue_ids unique, and exactly the protein block kept."""
    # 6x6 PAE with a recognizable pattern: pae[i, j] = 10*i + j.
    pae = np.fromfunction(lambda i, j: 10 * i + j, (6, 6), dtype=float)
    # chain A = protein (res 1,2,3); B/C/D = Mg (res 1 each, COLLIDES with protein 1); E = POP res 1.
    res_ids = [1, 2, 3, 1, 1, 1]
    chain_ids = ["A", "A", "A", "B", "C", "D"]
    with tempfile.TemporaryDirectory() as tmp:
        conf = _write_job(tmp, pae=pae, token_res_ids=res_ids, token_chain_ids=chain_ids)
        out = os.path.join(tmp, "job_pae.npz")
        L = confidences_to_npz(conf, out, job_id="job")
        assert L == 3, f"expected 3 protein tokens after restriction, got {L}"
        with np.load(out) as npz:
            got_pae = np.asarray(npz["pae"], dtype=float)
            got_ids = np.asarray(npz["residue_ids"]).astype(int)
        assert got_pae.shape == (3, 3), got_pae.shape
        assert list(got_ids) == [1, 2, 3], list(got_ids)
        assert np.unique(got_ids).shape[0] == got_ids.shape[0], "residue_ids must be unique"
        # The kept block must be exactly the protein-chain rows/cols (indices 0,1,2).
        assert np.allclose(got_pae, pae[:3, :3]), got_pae
    print("  ok  test_holo_multichain_restricted_to_protein_chain")


def test_single_chain_unchanged():
    """A single-chain fold is stored verbatim (no slicing)."""
    pae = np.fromfunction(lambda i, j: 10 * i + j, (4, 4), dtype=float)
    with tempfile.TemporaryDirectory() as tmp:
        conf = _write_job(
            tmp, pae=pae, token_res_ids=[1, 2, 3, 4], token_chain_ids=["A", "A", "A", "A"]
        )
        out = os.path.join(tmp, "job_pae.npz")
        L = confidences_to_npz(conf, out, job_id="job")
        assert L == 4, L
        with np.load(out) as npz:
            assert np.allclose(np.asarray(npz["pae"], dtype=float), pae)
            assert list(np.asarray(npz["residue_ids"]).astype(int)) == [1, 2, 3, 4]
            assert abs(float(npz["ptm"]) - 0.8) < 1e-6
            assert np.isnan(float(npz["iptm"]))  # iptm null -> NaN
    print("  ok  test_single_chain_unchanged")


def test_largest_chain_is_the_protein():
    """Protein chain need not be 'A' / first -- the chain with the most tokens wins."""
    pae = np.fromfunction(lambda i, j: 10 * i + j, (5, 5), dtype=float)
    # chain Z is the ion (1 token, first); chain A is the protein (4 tokens).
    res_ids = [1, 1, 2, 3, 4]
    chain_ids = ["Z", "A", "A", "A", "A"]
    with tempfile.TemporaryDirectory() as tmp:
        conf = _write_job(tmp, pae=pae, token_res_ids=res_ids, token_chain_ids=chain_ids)
        out = os.path.join(tmp, "job_pae.npz")
        L = confidences_to_npz(conf, out, job_id="job")
        assert L == 4, L
        with np.load(out) as npz:
            got_pae = np.asarray(npz["pae"], dtype=float)
            got_ids = np.asarray(npz["residue_ids"]).astype(int)
        assert list(got_ids) == [1, 2, 3, 4], list(got_ids)
        assert np.allclose(got_pae, pae[1:, 1:]), got_pae  # rows/cols 1..4 (chain A)
    print("  ok  test_largest_chain_is_the_protein")


def main() -> None:
    test_holo_multichain_restricted_to_protein_chain()
    test_single_chain_unchanged()
    test_largest_chain_is_the_protein()
    print("All 3 tests passed.")


if __name__ == "__main__":
    main()
