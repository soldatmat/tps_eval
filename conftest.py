"""Session-wide pytest guards for tps_eval's test suite.

Several tool tests fake out an external binary by assigning onto the module-level
``subprocess`` reference of the module under test, e.g.

    sws.subprocess.Popen = _FakePopen   # homology_search/test_swissprot_search.py

``sws.subprocess`` IS the one shared ``subprocess`` module, so that assignment is
process-wide and permanent — every later test in the same pytest session that
legitimately shells out silently gets the fake instead. That is why
``test_pocket_descriptors`` and ``selection::test_diversity_dedup_per_group`` failed
only when run as part of the full suite and passed in isolation (the pocket
traceback showed ``test_swissprot_search._FakePopen`` being handed an ``fpocket``
command).

Rather than rewriting five test files — which also have standalone ``main()`` shims
that must keep working outside pytest — snapshot and restore the globals around
every test. Tests keep patching however they like; the blast radius stops at the
test boundary.
"""

import subprocess

import pytest


@pytest.fixture(autouse=True)
def _restore_subprocess_globals():
    saved_run, saved_popen, saved_check = (
        subprocess.run,
        subprocess.Popen,
        subprocess.check_output,
    )
    try:
        yield
    finally:
        subprocess.run = saved_run
        subprocess.Popen = saved_popen
        subprocess.check_output = saved_check
