"""2D dimensionality-reduction backends for landscape maps.

Each function returns an (n, 2) coordinate array (PCA/PCoA also return the
%-variance of the two axes). Backends accept either a feature matrix X (n, d)
or, where noted, a precomputed (n, n) distance matrix — so the same renderer
serves embedding spaces (ESM/SaProt/EE) and similarity matrices (mmseqs /
foldseek).

t-SNE/UMAP/PaCMAP are imported lazily so the package imports without them.
"""
from __future__ import annotations

import numpy as np


def _zscore(X: np.ndarray) -> np.ndarray:
    mu = X.mean(0, keepdims=True)
    sd = X.std(0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def pca_2d(X: np.ndarray, zscore: bool = True):
    if zscore:
        X = _zscore(X)
    Xc = X - X.mean(0, keepdims=True)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    coords = U[:, :2] * S[:2]
    vr = S ** 2 / (S ** 2).sum()
    return coords, (100 * vr[0], 100 * vr[1])


def tsne_2d(data, precomputed: bool = False, perplexity: int = 30,
            random_state: int = 0, zscore: bool = True):
    from sklearn.manifold import TSNE
    if precomputed:
        return TSNE(n_components=2, metric="precomputed", init="random",
                    perplexity=perplexity, random_state=random_state).fit_transform(data)
    X = _zscore(data) if zscore else data
    Xc = X - X.mean(0, keepdims=True)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    n = min(50, Xc.shape[1])
    X50 = U[:, :n] * S[:n]
    return TSNE(n_components=2, init="pca", perplexity=perplexity,
                random_state=random_state).fit_transform(X50)


def umap_2d(data, precomputed: bool = False, n_neighbors: int = 15,
            min_dist: float = 0.1, random_state: int = 0, zscore: bool = True):
    import umap
    if precomputed:
        return umap.UMAP(n_components=2, metric="precomputed", n_neighbors=n_neighbors,
                         min_dist=min_dist, random_state=random_state).fit_transform(data)
    X = _zscore(data) if zscore else data
    return umap.UMAP(n_components=2, metric="euclidean", n_neighbors=n_neighbors,
                     min_dist=min_dist, random_state=random_state).fit_transform(X)


def pacmap_2d(X: np.ndarray, random_state: int = 0, zscore: bool = True):
    import pacmap
    if zscore:
        X = _zscore(X)
    return pacmap.PaCMAP(n_components=2, random_state=random_state).fit_transform(X)


def pcoa_2d(D: np.ndarray):
    """Classical MDS (principal coordinates) on a distance matrix D (n, n)."""
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    coords = V[:, :2] * np.sqrt(np.maximum(w[:2], 0.0))
    pos = w[w > 0].sum()
    pct = (100 * w[:2] / pos) if pos > 0 else np.array([0.0, 0.0])
    return coords, (pct[0], pct[1])
