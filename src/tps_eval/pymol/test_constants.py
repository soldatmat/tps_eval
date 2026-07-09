from __future__ import annotations

"""Self-contained tests for pymol/constants.py (PyMOL render settings).

Run from this directory:
    cd src/pymol && python test_constants.py
or under pytest:
    cd src/pymol && python -m pytest test_constants.py -q

Trivial sanity locks on the render constants: positive image width / DPI, a
0-or-1 ray flag (PyMOL's cmd.png ray is a boolean-ish int), and a positive
cartoon width. These guard against a sign/type slip that would break rendering.
"""

import tps_eval.pymol.constants as K


def test_png_width_positive_int():
    assert isinstance(K.PNG_WIDTH, int) and K.PNG_WIDTH > 0


def test_png_dpi_positive_int():
    assert isinstance(K.PNG_DPI, int) and K.PNG_DPI > 0


def test_png_ray_is_flag():
    assert K.PNG_RAY in (0, 1)


def test_secondary_structure_width_positive():
    assert isinstance(K.SMALL_SECONDARY_STRUCTURE_WIDTH, (int, float))
    assert K.SMALL_SECONDARY_STRUCTURE_WIDTH > 0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
