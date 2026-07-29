"""Canonical on-disk anchors for the src-layout repo — import these, don't recount `..`.

The package lives at ``<repo>/src/tps_eval/``, so a module at
``<repo>/src/tps_eval/<subdir>/<module>.py`` is FOUR levels below the repo root.
Hand-rolled ``Path(__file__).parent.parent.parent`` chains got this wrong after the
src-layout migration and pointed at ``<repo>/src`` instead — which silently broke
every ``vendor/`` shell-out (``proteinmpnn_score`` and ``self_consistency`` wrote an
all-NaN column for two production runs; ``catapro``'s predictor and the dashboard's
default reference bands were dead the same way, and the AF3 fan-out looked for
``submit_job.sh`` under ``src/scripts/``).

Only ``tps_eval`` package data (reference-stats JSONs, the accession list, HTML
templates) is addressed from `PACKAGE_DIR`; anything OUTSIDE the package — `vendor/`
submodules, `scripts/`, the gitignored `data/` — hangs off `REPO_ROOT`.

This assumes the package is used from a checkout (editable install or the
``PYTHONPATH=<repo>/src`` the tool wrappers export), which is how every cluster runs
it. A non-editable wheel copied into site-packages has no repo around it and no
``vendor/`` to find either way.
"""

from pathlib import Path

#: ``<repo>/src/tps_eval`` — the import package itself (package data lives here).
PACKAGE_DIR = Path(__file__).resolve().parent
#: ``<repo>/src`` — the src-layout source root.
SRC_DIR = PACKAGE_DIR.parent
#: ``<repo>`` — the checkout root.
REPO_ROOT = SRC_DIR.parent
#: ``<repo>/vendor`` — git submodules (ProteinMPNN, CataPro, TmProt, aggrescan3d, …).
VENDOR_DIR = REPO_ROOT / "vendor"
#: ``<repo>/scripts`` — wrappers, job scripts, orchestrators.
SCRIPTS_DIR = REPO_ROOT / "scripts"
#: ``<repo>/src/tps_eval/reference_stats`` — committed MARTS-DB reference bands.
REFERENCE_STATS_DIR = PACKAGE_DIR / "reference_stats"
