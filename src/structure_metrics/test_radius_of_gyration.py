from __future__ import annotations

"""Self-contained tests for radius_of_gyration.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_radius_of_gyration.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_radius_of_gyration.py -q

The geometric tests use synthetic point sets whose Rg / gyration-tensor
eigenvalues are known in closed form, so the assertions are exact (no fold model
needed). The I/O tests write tiny PDB files to a temp dir and exercise the
flat-directory collection, ID keying, and NaN-on-failure contract.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from radius_of_gyration import (
    COLUMNS,
    ca_coordinates,
    gyration_metrics,
    radius_of_gyration_dir,
)


def _approx(a, b, tol=1e-6):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_empty_is_nan():
    m = gyration_metrics(np.empty((0, 3)))
    assert m["n_residues"] == 0
    for k in ("radius_of_gyration", "asphericity", "acylindricity",
              "principal_radius_1", "principal_radius_2", "principal_radius_3"):
        assert np.isnan(m[k]), k


def test_single_point_zero():
    m = gyration_metrics(np.zeros((1, 3)))
    assert m["n_residues"] == 1
    _approx(m["radius_of_gyration"], 0.0)
    _approx(m["asphericity"], 0.0)
    _approx(m["acylindricity"], 0.0)


def test_points_on_x_axis():
    # Two points at +-d along x: com at origin, Rg = d.
    d = 5.0
    coords = np.array([[-d, 0, 0], [d, 0, 0]], dtype=float)
    m = gyration_metrics(coords)
    _approx(m["radius_of_gyration"], d)
    # Gyration tensor: λ1 = d^2 (along x), λ2 = λ3 = 0.
    _approx(m["principal_radius_1"], d)
    _approx(m["principal_radius_2"], 0.0)
    _approx(m["principal_radius_3"], 0.0)
    # asphericity = λ1 - (λ2+λ3)/2 = d^2; acylindricity = λ2 - λ3 = 0.
    _approx(m["asphericity"], d * d)
    _approx(m["acylindricity"], 0.0)
    # Rg^2 == λ1+λ2+λ3.
    rg2 = m["radius_of_gyration"] ** 2
    lam_sum = (m["principal_radius_1"] ** 2 + m["principal_radius_2"] ** 2
               + m["principal_radius_3"] ** 2)
    _approx(rg2, lam_sum)


def test_isotropic_sphere_low_asphericity():
    # 6 points on the +-axes at radius r: perfectly isotropic gyration tensor.
    r = 3.0
    coords = np.array([
        [r, 0, 0], [-r, 0, 0],
        [0, r, 0], [0, -r, 0],
        [0, 0, r], [0, 0, -r],
    ], dtype=float)
    m = gyration_metrics(coords)
    # Each eigenvalue = (2 r^2)/6 = r^2/3; Rg^2 = r^2.
    _approx(m["radius_of_gyration"], r)
    _approx(m["asphericity"], 0.0)
    _approx(m["acylindricity"], 0.0)
    _approx(m["principal_radius_1"], r / np.sqrt(3.0))


def test_translation_invariance():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(50, 3)) * 4.0
    base = gyration_metrics(coords)
    shifted = gyration_metrics(coords + np.array([100.0, -50.0, 7.0]))
    for k in COLUMNS[1:-1]:
        _approx(base[k], shifted[k], tol=1e-6)
    # Eigenvalue ordering: λ1 >= λ2 >= λ3.
    assert base["principal_radius_1"] >= base["principal_radius_2"] >= base["principal_radius_3"]
    assert base["asphericity"] >= 0.0
    assert base["acylindricity"] >= -1e-9


def _write_pdb(path, coords, chain="A"):
    """Write a minimal PDB with one CA per residue (ALA) at the given coords."""
    with open(path, "w") as fh:
        for i, (x, y, z) in enumerate(coords, start=1):
            fh.write(
                f"ATOM  {i:>5d}  CA  ALA {chain}{i:>4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        fh.write("END\n")


def test_ca_coordinates_and_hetatm_skipped():
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(20, 3)) * 6.0
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "design.pdb")
        with open(p, "w") as fh:
            for i, (x, y, z) in enumerate(coords, start=1):
                fh.write(
                    f"ATOM  {i:>5d}  CA  ALA A{i:>4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
                )
            # A HETATM calcium ion whose atom is named CA must be ignored.
            fh.write(
                "HETATM  999 CA    CA A 999     999.000 999.000 999.000  1.00  0.00          CA\n"
            )
            fh.write("END\n")
        got = ca_coordinates(p)
        assert got.shape == (20, 3), got.shape
        np.testing.assert_allclose(np.sort(got, axis=0), np.sort(coords, axis=0), atol=1e-3)


def test_dir_keyed_by_id_and_nan_on_failure():
    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        os.makedirs(structs)
        # Good structure: two CA along x, d=5 -> Rg=5.
        _write_pdb(os.path.join(structs, "good.pdb"),
                   np.array([[-5, 0, 0], [5, 0, 0]], dtype=float))
        # Broken structure: no parseable ATOM records -> NaN row, must not abort.
        with open(os.path.join(structs, "broken.pdb"), "w") as fh:
            fh.write("this is not a pdb\n")
        df = radius_of_gyration_dir(structs)

        assert list(df.columns) == COLUMNS
        assert set(df["ID"]) == {"good", "broken"}
        # Default save path is the sibling CSV.
        assert os.path.isfile(structs + "_radius_of_gyration.csv")

        good = df.set_index("ID").loc["good"]
        _approx(float(good["radius_of_gyration"]), 5.0, tol=1e-3)
        assert int(good["n_residues"]) == 2

        broken = df.set_index("ID").loc["broken"]
        assert np.isnan(broken["radius_of_gyration"])
        assert int(broken["n_residues"]) == 0


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
