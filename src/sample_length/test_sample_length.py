from __future__ import annotations

"""Self-contained tests for sample_length.py (length-distribution sampler).

Run from this directory:
    cd src/sample_length && python test_sample_length.py
or under pytest:
    cd src/sample_length && python -m pytest test_sample_length.py -q

Writes a tiny synthetic FASTA to a temp dir, builds the empirical length sampler, and
asserts: sampled lengths are drawn only from the observed lengths, the sampler is
reproducible under a fixed np.random.seed, the empirical frequencies converge to the
input length proportions on a large draw, and the return_counts=True variant returns
(length, count) pairs summing to n.
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sample_length import get_sampler  # noqa: E402


def _write_fasta(path, lengths):
    """One record per entry, sequence of the given length (poly-A)."""
    with open(path, "w") as fh:
        for i, L in enumerate(lengths):
            fh.write(f">seq{i}\n{'A' * L}\n")


def test_sampled_lengths_within_support():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [10, 10, 20, 30])
        sampler = get_sampler(fa)
        np.random.seed(0)
        out = sampler(100)
        assert len(out) == 100
        assert set(np.unique(out)).issubset({10, 20, 30})


def test_reproducible_with_seed():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [5, 7, 7, 11, 11, 11])
        sampler = get_sampler(fa)
        np.random.seed(42)
        a = sampler(50)
        np.random.seed(42)
        b = sampler(50)
        assert list(a) == list(b)


def test_frequencies_track_input_proportions():
    """A length appearing 3x as often should be sampled ~3x as often (large draw)."""
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        # length 100 appears 3x, length 200 appears 1x -> ~75% / ~25%.
        _write_fasta(fa, [100, 100, 100, 200])
        sampler = get_sampler(fa)
        np.random.seed(1)
        out = sampler(20000)
        frac_100 = np.mean(out == 100)
        assert abs(frac_100 - 0.75) < 0.03, frac_100


def test_return_counts_variant():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "in.fasta")
        _write_fasta(fa, [10, 20, 20, 30])
        sampler = get_sampler(fa, return_counts=True)
        np.random.seed(7)
        pairs = sampler(200)
        # list of (length, count) pairs
        assert all(len(p) == 2 for p in pairs)
        lengths = {int(l) for l, _ in pairs}
        assert lengths.issubset({10, 20, 30})
        assert sum(int(c) for _, c in pairs) == 200


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
