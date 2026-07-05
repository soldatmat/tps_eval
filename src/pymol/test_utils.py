from __future__ import annotations

"""Self-contained tests for pymol/utils.py.

Run from this directory:
    cd src/pymol && python test_utils.py
or under pytest:
    cd src/pymol && python -m pytest test_utils.py -q

`show_organic_and_metals` takes the PyMOL `cmd` module as an argument and only
issues show/color calls, so it is testable WITHOUT PyMOL by passing a recording
stub. We assert the exact sequence of PyMOL commands (selections + representations)
it emits — the part that would silently break a render if reordered/mis-selected.
"""


class _CmdRecorder:
    """Minimal stand-in for pymol.cmd that records show()/color() calls."""

    def __init__(self):
        self.calls = []

    def show(self, representation, selection):
        self.calls.append(("show", representation, selection))

    def color(self, color, selection):
        self.calls.append(("color", color, selection))


from utils import show_organic_and_metals


def test_show_organic_and_metals_emits_expected_calls():
    rec = _CmdRecorder()
    show_organic_and_metals(rec)
    assert rec.calls == [
        ("show", "sticks", "structure and organic"),
        ("color", "atomic", "structure and organic"),
        ("color", "gray", "structure and organic and elem C"),
        ("show", "spheres", "structure and metals"),
    ]


def test_carbons_recolored_after_atomic():
    """The gray-carbon override must come AFTER the atomic coloring, else it is
    clobbered (this ordering is the whole point of the helper)."""
    rec = _CmdRecorder()
    show_organic_and_metals(rec)
    atomic_i = rec.calls.index(("color", "atomic", "structure and organic"))
    gray_i = rec.calls.index(("color", "gray", "structure and organic and elem C"))
    assert gray_i > atomic_i


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
