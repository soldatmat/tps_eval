from __future__ import annotations

"""Self-contained tests for the PURE helpers in esmfold.py.

Run from this directory:
    cd src/esmfold && python test_esmfold.py
or under pytest:
    cd src/esmfold && python -m pytest test_esmfold.py -q

ESMFold itself (torch/transformers) is NOT exercised here (NEEDS-AURUM/GPU). We test
only the pure, model-free helpers, which is where the correctness-critical contracts
live: the B-factor pLDDT 0-1 -> 0-100 rescale (the drop-in-with-AlphaFold invariant),
the PDB CA residue-id axis reader, the id sanitizer, and the sibling PAE-dir path. The
PAE/pTM extractors are exercised with a tiny numpy-backed fake tensor (no torch), which
validates the batch-slice + crop-to-length logic and the graceful None handling.
"""

import os
import sys

import numpy as np

# esmfold/ is not a package; import the flat module by putting its own dir first
# (mirrors test_extract_pae.py). esmfold.py adds src/ to sys.path itself for its
# `from data.sequences import ...` at import time.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from esmfold import (  # noqa: E402
    _default_pae_dir,
    _extract_pae,
    _extract_ptm,
    _rescale_bfactor_to_0_100,
    _residue_ids_from_pdb,
    _sanitize_id,
)


class _FakeTensor:
    """Mimics the tiny slice of the torch tensor API the extractors touch."""

    def __init__(self, arr):
        self._arr = np.asarray(arr)

    def detach(self):
        return self

    def to(self, *a, **k):
        return self

    def float(self):
        return self

    def numpy(self):
        return self._arr


def _atom_line(serial, atom, resname, chain, resseq, x, y, z, bfac, element, icode=" "):
    return (
        f"ATOM  {serial:>5d} {atom:<4s}{resname:>3s} {chain}{resseq:>4d}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfac:6.2f}          {element:>2s}"
    )


def test_rescale_bfactor_0_1_to_0_100():
    """ESMFold writes pLDDT in [0,1]; the rescale must multiply the B-factor by 100."""
    line = _atom_line(1, " CA ", "ALA", "A", 1, 1.0, 2.0, 3.0, 0.85, "C")
    out = _rescale_bfactor_to_0_100(line + "\nEND")
    b = float(out.splitlines()[0][60:66])
    assert abs(b - 85.0) < 1e-6, b
    # column layout preserved: coords untouched.
    assert out.splitlines()[0][:60] == line[:60]


def test_rescale_handles_hetatm_and_full_range():
    het = "HETATM" + _atom_line(2, " MG ", " MG", "B", 1, 0.0, 0.0, 0.0, 1.00, "MG")[6:]
    out = _rescale_bfactor_to_0_100(het)
    assert abs(float(out.splitlines()[0][60:66]) - 100.0) < 1e-6


def test_rescale_leaves_non_atom_lines_untouched():
    txt = "REMARK something\nTER\nEND\n"
    assert _rescale_bfactor_to_0_100(txt).startswith("REMARK something")


def test_residue_ids_from_pdb_reads_author_numbering():
    lines = [
        _atom_line(1, " N  ", "ALA", "A", 5, 0, 0, 0, 50.0, "N"),
        _atom_line(2, " CA ", "ALA", "A", 5, 0, 0, 0, 50.0, "C"),
        _atom_line(3, " CA ", "GLY", "A", 6, 0, 0, 0, 60.0, "C"),
        _atom_line(4, " CA ", "SER", "A", 8, 0, 0, 0, 70.0, "C"),  # gap 7 (non-contiguous)
    ]
    ids = _residue_ids_from_pdb("\n".join(lines) + "\n")
    assert ids == [5, 6, 8], ids  # one entry per CA, author numbers, gaps preserved


def test_sanitize_id():
    assert _sanitize_id("design_42 description here") == "design_42"
    assert _sanitize_id("  spaced  ") == "spaced"
    assert _sanitize_id("a/b/c") == "a_b_c"


def test_default_pae_dir_is_sibling():
    assert _default_pae_dir("/x/y/gen_esmfold_structs") == "/x/y/gen_esmfold_structs_pae"
    assert _default_pae_dir("/x/y/structs/") == "/x/y/structs_pae"


def test_extract_pae_batch_slice_and_crop():
    # (batch, L, L) with L=4; seq_len=3 -> crop to 3x3, batch element 0 taken.
    arr = np.arange(2 * 4 * 4, dtype=np.float32).reshape(2, 4, 4)
    out = _extract_pae({"predicted_aligned_error": _FakeTensor(arr)}, seq_len=3)
    assert out.shape == (3, 3), out.shape
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, arr[0, :3, :3])


def test_extract_pae_missing_returns_none():
    assert _extract_pae({}, seq_len=5) is None


def test_extract_ptm_scalar():
    assert abs(_extract_ptm({"ptm": _FakeTensor(np.array(0.73, dtype=np.float32))}) - 0.73) < 1e-6
    assert _extract_ptm({}) is None


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
