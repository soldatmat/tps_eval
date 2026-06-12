"""tps_eval shared data-loading helpers (sequences, embeddings, results).

This package intentionally carries an ``__init__.py`` (rather than relying on
implicit namespace-package resolution): some satellite conda envs that reuse
this repo — e.g. an ESMFold/ProteinMPNN env that also ships its own top-level
``data`` package — would otherwise shadow ``src/data`` (a regular package on a
later sys.path entry beats a namespace portion earlier on the path). Making this
a regular package guarantees ``from data.sequences import ...`` resolves here
once ``src/`` is on sys.path, regardless of the active env.
"""
