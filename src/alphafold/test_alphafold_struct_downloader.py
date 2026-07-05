from __future__ import annotations

"""Self-contained tests for alphafold_struct_downloader.py (AFDB URL / save-path
construction + graceful failure).

Run from this directory:
    cd src/alphafold && python test_alphafold_struct_downloader.py
or under pytest:
    cd src/alphafold && python -m pytest test_alphafold_struct_downloader.py -q

NO network: we monkeypatch the module's `requests.get` with a stub that records the
requested URL and returns a fake response, then assert on the constructed AFDB v6 URL,
the (accession vs save-name) file naming, the 404 (not-found) no-write path, and that
a raised request no longer crashes the finally block (regression: `response` was
referenced unbound). Uses tempfile only.
"""

import os
import sys
import tempfile
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import alphafold.alphafold_struct_downloader as dl


class _FakeResponse:
    def __init__(self, status_code=200, content=b"PDBDATA"):
        self.status_code = status_code
        self.content = content


def _install_fake_get(recorded, response=None, raise_exc=None):
    orig = dl.requests.get

    def fake_get(url, *a, **k):
        recorded.append(url)
        if raise_exc is not None:
            raise raise_exc
        return response if response is not None else _FakeResponse()

    dl.requests.get = fake_get
    return lambda: setattr(dl.requests, "get", orig)


def test_url_construction_and_write():
    recorded = []
    restore = _install_fake_get(recorded, _FakeResponse(200, b"COORDS"))
    try:
        with tempfile.TemporaryDirectory() as d:
            dl.download_af_struct("P12345", Path(d))
            assert recorded == [
                "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.pdb"
            ], recorded
            out = Path(d) / "P12345.pdb"
            assert out.is_file()
            assert out.read_bytes() == b"COORDS"
    finally:
        restore()


def test_tuple_saves_under_save_name():
    """(uniprot_id, save_name): URL uses the accession, file uses the save name."""
    recorded = []
    restore = _install_fake_get(recorded, _FakeResponse(200, b"X"))
    try:
        with tempfile.TemporaryDirectory() as d:
            dl.download_af_struct(("Q99999", "my_enzyme"), Path(d))
            assert "AF-Q99999-F1-model_v6.pdb" in recorded[0]
            assert (Path(d) / "my_enzyme.pdb").is_file()
            assert not (Path(d) / "Q99999.pdb").exists()
    finally:
        restore()


def test_not_found_writes_nothing():
    recorded = []
    restore = _install_fake_get(recorded, _FakeResponse(404, b""))
    try:
        with tempfile.TemporaryDirectory() as d:
            dl.download_af_struct("NOPE", Path(d))
            assert list(Path(d).iterdir()) == []
    finally:
        restore()


def test_request_exception_does_not_crash():
    """A raised request must not blow up the finally with an UnboundLocalError.
    (Regression: `response` was referenced in finally even when get() raised.)"""
    recorded = []
    restore = _install_fake_get(recorded, raise_exc=ConnectionError("boom"))
    try:
        with tempfile.TemporaryDirectory() as d:
            # max_fails_count=0 so it does not recurse/retry (which would re-raise).
            dl.download_af_struct("P00000", Path(d), fails_count=99, max_fails_count=3)
            assert list(Path(d).iterdir()) == []
    finally:
        restore()


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
