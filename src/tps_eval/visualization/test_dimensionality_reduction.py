from __future__ import annotations

"""Self-contained tests for visualization/dimensionality_reduction.py.

Run from this directory (so the package-relative import resolves like the
package does):
    cd src && python -m visualization.test_dimensionality_reduction
or under pytest:
    cd src && python -m pytest visualization/test_dimensionality_reduction.py -q

Only the numpy-only backends (PCA, PCoA, z-score) are exercised — they are pure
and deterministic. t-SNE / UMAP / PaCMAP are lazy-imported wrappers over external
libraries and are not tested here (no rendering, no stochastic backends).
Assertions use synthetic point sets whose geometry is known in closed form.
"""

import sys
from pathlib import Path


import numpy as np

from tps_eval.visualization.dimensionality_reduction import _zscore, pca_2d, pcoa_2d


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_zscore_zero_mean_unit_std():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 4)) * 3.0 + 5.0
    Z = _zscore(X)
    np.testing.assert_allclose(Z.mean(0), 0.0, atol=1e-9)
    np.testing.assert_allclose(Z.std(0), 1.0, atol=1e-9)


def test_zscore_constant_column_no_divide_by_zero():
    # A constant column has sd 0; the guard sets sd->1 so output is finite (all 0).
    X = np.array([[1.0, 7.0], [2.0, 7.0], [3.0, 7.0]])
    Z = _zscore(X)
    assert np.all(np.isfinite(Z))
    np.testing.assert_allclose(Z[:, 1], 0.0)


def test_pca_2d_shape_and_variance():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(50, 6))
    coords, (v0, v1) = pca_2d(X, zscore=True)
    assert coords.shape == (50, 2)
    # Variance percentages are ordered and in (0, 100].
    assert 0.0 < v1 <= v0 <= 100.0


def test_pca_2d_recovers_dominant_axis():
    # Data spread mostly along axis 0 -> PC1 should capture most variance.
    rng = np.random.default_rng(2)
    n = 300
    X = np.column_stack([
        rng.normal(scale=10.0, size=n),   # big spread
        rng.normal(scale=0.1, size=n),    # tiny spread
    ])
    _, (v0, v1) = pca_2d(X, zscore=False)
    assert v0 > v1
    assert v0 > 90.0  # first component dominates


def test_pcoa_2d_preserves_pairwise_distances():
    # PCoA on the Euclidean distance matrix of a planar point set recovers the
    # configuration up to rotation/reflection -> pairwise distances preserved.
    pts = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 4.0], [3.0, 4.0], [1.5, 2.0]])
    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    coords, (p0, p1) = pcoa_2d(D)
    assert coords.shape == (5, 2)
    D_rec = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.testing.assert_allclose(np.sort(D_rec.ravel()), np.sort(D.ravel()), atol=1e-6)
    # Two axes of a genuinely 2D set explain ~all positive variance.
    _approx(p0 + p1, 100.0, tol=1e-4)


def test_pcoa_2d_variance_percentages_ordered():
    rng = np.random.default_rng(3)
    pts = rng.normal(size=(20, 5))
    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    _, (p0, p1) = pcoa_2d(D)
    assert p0 >= p1 >= 0.0


def test_pca_2d_single_feature_pads_second_axis():
    # A 1-feature matrix has only one principal component: coords must still be (n, 2)
    # (2nd axis padded with zeros) and the 2nd variance % is NaN (undefined) -- not an
    # IndexError on vr[1] (regression).
    X = np.linspace(0.0, 1.0, 20).reshape(-1, 1)
    coords, (v0, v1) = pca_2d(X, zscore=False)
    assert coords.shape == (20, 2)
    np.testing.assert_allclose(coords[:, 1], 0.0)
    _approx(v0, 100.0)
    assert np.isnan(v1)


def test_pca_2d_degenerate_input_no_nan_variance():
    # All-identical rows -> zero spread. The variance % must degrade to 0.0 (guarded),
    # not NaN from a 0/0 division.
    X = np.ones((10, 4))
    _, (v0, v1) = pca_2d(X, zscore=False)
    assert v0 == 0.0 and v1 == 0.0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
