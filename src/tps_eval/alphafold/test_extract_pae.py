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
from pathlib import Path

import numpy as np

from tps_eval.alphafold.extract_pae import (  # noqa: E402
    _coerce_scalar,
    _load_ptm_iptm,
    _resolve_af_output,
    confidences_to_npz,
    extract_af_output,
)


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


def test_coerce_scalar():
    """ptm/iptm coercion: null -> NaN, list -> mean (ignoring nulls), garbage -> NaN."""
    assert np.isnan(_coerce_scalar(None))
    assert abs(_coerce_scalar(0.8) - 0.8) < 1e-6
    assert abs(_coerce_scalar([1.0, 3.0]) - 2.0) < 1e-6
    assert abs(_coerce_scalar([1.0, None, 3.0]) - 2.0) < 1e-6  # per-chain list w/ a null
    assert np.isnan(_coerce_scalar([]))
    assert np.isnan(_coerce_scalar("not-a-number"))
    print("  ok  test_coerce_scalar")


def test_res_ids_absent_defaults_to_arange():
    """When token_res_ids is missing, residue_ids falls back to 1..L."""
    pae = np.zeros((3, 3), dtype=float)
    with tempfile.TemporaryDirectory() as tmp:
        conf = os.path.join(tmp, "job_confidences.json")
        with open(conf, "w") as fh:
            json.dump({"pae": pae.tolist()}, fh)  # no token_res_ids / chain_ids
        out = os.path.join(tmp, "job_pae.npz")
        L = confidences_to_npz(conf, out, job_id="job")
        assert L == 3
        with np.load(out) as npz:
            assert list(np.asarray(npz["residue_ids"]).astype(int)) == [1, 2, 3]
    print("  ok  test_res_ids_absent_defaults_to_arange")


def test_non_square_and_missing_pae_raise():
    with tempfile.TemporaryDirectory() as tmp:
        # non-square PAE -> ValueError
        conf = os.path.join(tmp, "a_confidences.json")
        with open(conf, "w") as fh:
            json.dump({"pae": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]}, fh)
        try:
            confidences_to_npz(conf, os.path.join(tmp, "a.npz"), job_id="a")
            raise AssertionError("expected ValueError for non-square PAE")
        except ValueError:
            pass
        # missing 'pae' field -> KeyError
        conf2 = os.path.join(tmp, "b_confidences.json")
        with open(conf2, "w") as fh:
            json.dump({"token_res_ids": [1, 2]}, fh)
        try:
            confidences_to_npz(conf2, os.path.join(tmp, "b.npz"), job_id="b")
            raise AssertionError("expected KeyError for missing pae")
        except KeyError:
            pass
    print("  ok  test_non_square_and_missing_pae_raise")


def test_res_ids_length_mismatch_raises():
    with tempfile.TemporaryDirectory() as tmp:
        conf = os.path.join(tmp, "job_confidences.json")
        with open(conf, "w") as fh:
            json.dump({"pae": np.zeros((3, 3)).tolist(), "token_res_ids": [1, 2]}, fh)
        try:
            confidences_to_npz(conf, os.path.join(tmp, "j.npz"), job_id="job")
            raise AssertionError("expected ValueError for res_ids/PAE length mismatch")
        except ValueError:
            pass
    print("  ok  test_res_ids_length_mismatch_raises")


def test_load_ptm_iptm_falls_back_to_matrix_file():
    """When no *_summary_confidences.json exists, the scalars are read from the big
    matrix file itself if a given AF3 version inlined them there."""
    with tempfile.TemporaryDirectory() as tmp:
        conf = os.path.join(tmp, "job_confidences.json")
        with open(conf, "w") as fh:
            json.dump({"pae": [[0.0]], "ptm": 0.9, "iptm": None}, fh)
        ptm, iptm = _load_ptm_iptm(conf, "job")
        assert abs(ptm - 0.9) < 1e-6
        assert np.isnan(iptm)
    print("  ok  test_load_ptm_iptm_falls_back_to_matrix_file")


def test_resolve_af_output_from_structs_sibling():
    """--structs_dir with a sibling af_output/ holding <job>/<job>_confidences.json."""
    with tempfile.TemporaryDirectory() as tmp:
        af_output = os.path.join(tmp, "af_output")
        job_dir = os.path.join(af_output, "d1")
        os.makedirs(job_dir)
        with open(os.path.join(job_dir, "d1_confidences.json"), "w") as fh:
            json.dump({"pae": [[0.0]]}, fh)
        structs = os.path.join(tmp, "structs")
        os.makedirs(structs)
        resolved = _resolve_af_output(None, structs)
        assert os.path.realpath(str(resolved)) == os.path.realpath(af_output)
    print("  ok  test_resolve_af_output_from_structs_sibling")


def test_extract_af_output_skips_existing():
    """extract_af_output writes one npz per job and counts a pre-existing npz as kept
    (skip_existing) without re-reading it."""
    with tempfile.TemporaryDirectory() as tmp:
        af_output = os.path.join(tmp, "af_output")
        for job in ("d1", "d2"):
            jd = os.path.join(af_output, job)
            os.makedirs(jd)
            with open(os.path.join(jd, f"{job}_confidences.json"), "w") as fh:
                json.dump({"pae": np.zeros((2, 2)).tolist()}, fh)
        pae_dir = os.path.join(tmp, "pae")
        n1 = extract_af_output(Path(af_output), pae_dir, skip_existing=True)
        assert n1 == 2
        assert os.path.isfile(os.path.join(pae_dir, "d1_pae.npz"))
        assert os.path.isfile(os.path.join(pae_dir, "d2_pae.npz"))
        # Second pass: both skipped but still counted.
        n2 = extract_af_output(Path(af_output), pae_dir, skip_existing=True)
        assert n2 == 2
    print("  ok  test_extract_af_output_skips_existing")


def main() -> None:
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
