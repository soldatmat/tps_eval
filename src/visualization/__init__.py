"""tps_eval.visualization — dataset-level landscape maps of a protein set.

Project a MARTS-DB-style protein set into 2D from any representation (PLM /
structure-aware / similarity-matrix / domain-feature) and render it coloured by
first-cyclization class, with the carbon-ordered substrate-type palette.

This is a *standalone* visualization layer (not wired into the per-design
evaluation orchestrator): representations are produced by the per-domain
packages (src/esm, src/saprot, src/foldseek, src/enzyme_explorer, src/specificity,
…) and exported as feature-matrix / distance CSVs; this package consumes those
and lays them out (PCA / t-SNE / UMAP / PaCMAP / PCoA).

Distinct from:
  - src/plot/  — per-design metric-comparison charts (boxplot/density/categorical)
  - src/pymol/ — 3D structure renders
"""
from . import palette, dimensionality_reduction, landscape_map  # noqa: F401
