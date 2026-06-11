from __future__ import annotations

"""Self-contained tests for interdomain_pae.py.

Run from this directory (so the flat-module imports resolve like the runner does):
    cd src/structure_metrics && python test_interdomain_pae.py
or under pytest:
    cd src/structure_metrics && python -m pytest test_interdomain_pae.py -q

These tests use synthetic PAE matrices and synthetic EE region dicts (the
``residues_mapping``-keyed shape EE's JSON sidecar produces) so the inter-domain
block reduction is checked in closed form — no fold model and no EnzymeExplorer
needed. The directory driver is exercised against a tiny structs dir + pae dir.

NOTE: importing this module pulls in interdomain_pae, which imports
enzyme_explorer.domain_composition (for _structure_ids / detect_domains_json reuse).
That import does NOT require EnzymeExplorer at import time (EE is imported lazily
inside detect_domains_json), so these tests run in any env with numpy+pandas.
"""

import os
import tempfile

import numpy as np
import pandas as pd

from interdomain_pae import (
    BASE_COLUMNS,
    design_row,
    domain_residue_numbers,
    extract_interdomain_pae_dir,
    interdomain_pae_blocks,
    load_pae,
)


def _approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _region(module_id, domain, resis):
    """An EE-JSON-shaped region dict: residues_mapping keyed by str(resi)."""
    return {
        "module_id": module_id,
        "domain": domain,
        "tmscore": 0.9,
        "residues_mapping": {str(r): r for r in resis},
    }


def test_domain_residue_numbers_keys_by_module_id():
    regions = [
        _region("d_alpha_0", "alpha", [1, 2, 3]),
        _region("d_alpha_1", "alpha", [10, 11]),  # two alphas stay distinct
        {"module_id": "d_beta_0", "domain": "beta", "residues_mapping": {}},  # empty dropped
    ]
    out = domain_residue_numbers(regions)
    assert set(out) == {"d_alpha_0", "d_alpha_1"}
    assert out["d_alpha_0"] == [1, 2, 3]
    assert out["d_alpha_1"] == [10, 11]


def test_interdomain_block_two_domains_known_value():
    # PAE: low within a domain, high across. Residues 1-2 = domain A, 3-4 = domain B.
    # Build an asymmetric matrix so both-direction averaging is actually exercised.
    pae = np.array(
        [
            [0.5, 0.5, 20.0, 22.0],
            [0.5, 0.5, 18.0, 24.0],
            [10.0, 12.0, 0.5, 0.5],
            [14.0, 16.0, 0.5, 0.5],
        ],
        dtype=float,
    )
    res_ids = np.array([1, 2, 3, 4])
    domains = {"A": [1, 2], "B": [3, 4]}
    pairs = interdomain_pae_blocks(pae, res_ids, domains)
    assert set(pairs) == {"A_B"}
    mean_ab = pae[np.ix_([0, 1], [2, 3])].mean()  # 21.0
    mean_ba = pae[np.ix_([2, 3], [0, 1])].mean()  # 13.0
    _approx(pairs["A_B"], 0.5 * (mean_ab + mean_ba))  # 17.0
    _approx(pairs["A_B"], 17.0)


def test_three_domains_mean_and_max():
    pae = np.full((6, 6), 5.0)
    # boost the C-pair blocks so max picks them
    pae[4:6, 0:2] = 30.0
    pae[0:2, 4:6] = 30.0
    res_ids = np.array([1, 2, 3, 4, 5, 6])
    domains = {"A": [1, 2], "B": [3, 4], "C": [5, 6]}
    pairs = interdomain_pae_blocks(pae, res_ids, domains)
    assert set(pairs) == {"A_B", "A_C", "B_C"}
    _approx(pairs["A_B"], 5.0)
    _approx(pairs["A_C"], 30.0)  # both directions 30
    _approx(pairs["B_C"], 5.0)


def test_residue_ids_nontrivial_numbering():
    # PAE rows correspond to PDB resis 101,102,103,104 (not 0..3): mapping must use
    # the residue_ids axis, not the raw resi number.
    pae = np.array(
        [
            [0.0, 0.0, 9.0, 9.0],
            [0.0, 0.0, 9.0, 9.0],
            [9.0, 9.0, 0.0, 0.0],
            [9.0, 9.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    res_ids = np.array([101, 102, 103, 104])
    domains = {"A": [101, 102], "B": [103, 104]}
    pairs = interdomain_pae_blocks(pae, res_ids, domains)
    _approx(pairs["A_B"], 9.0)


def _write_npz(path, pae, res_ids):
    np.savez_compressed(
        path,
        pae=np.asarray(pae, dtype=np.float32),
        residue_ids=np.asarray(res_ids, dtype=np.int32),
        n_residues=np.int64(len(res_ids)),
        source="esmfold",
    )


def test_load_pae_missing_returns_none():
    pae, ids = load_pae("/no/such/file_pae.npz")
    assert pae is None and ids is None


def test_load_pae_default_axis_when_no_residue_ids():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x_pae.npz")
        np.savez_compressed(p, pae=np.zeros((3, 3), dtype=np.float32))
        pae, ids = load_pae(p)
        assert pae.shape == (3, 3)
        assert list(ids) == [1, 2, 3]


def test_design_row_single_domain_is_na():
    with tempfile.TemporaryDirectory() as d:
        _write_npz(os.path.join(d, "s_pae.npz"), np.ones((4, 4)), [1, 2, 3, 4])
        regions = [_region("s_alpha_0", "alpha", [1, 2, 3, 4])]
        row = design_row("s", d, regions)
        assert np.isnan(row["mean_interdomain_pae"])
        assert np.isnan(row["max_interdomain_pae"])
        assert row["n_domains"] == 1  # recorded even for N/A


def test_design_row_missing_pae_is_na_but_records_ndomains():
    with tempfile.TemporaryDirectory() as d:
        regions = [_region("s_alpha_0", "alpha", [1, 2]), _region("s_beta_0", "beta", [3, 4])]
        row = design_row("nope", d, regions)  # no npz on disk
        assert np.isnan(row["mean_interdomain_pae"])
        assert np.isnan(row["max_interdomain_pae"])
        assert row["n_domains"] == 2  # EE ran -> n_domains recorded even without PAE


def test_design_row_zero_domains_no_pae():
    with tempfile.TemporaryDirectory() as d:
        row = design_row("nope", d, [])  # no domains, no npz
        assert np.isnan(row["mean_interdomain_pae"])
        assert row["n_domains"] == 0


def test_design_row_two_domains_scored_and_per_pair():
    with tempfile.TemporaryDirectory() as d:
        pae = np.full((4, 4), 3.0)
        pae[0:2, 2:4] = 12.0
        pae[2:4, 0:2] = 8.0
        _write_npz(os.path.join(d, "t_pae.npz"), pae, [1, 2, 3, 4])
        regions = [_region("t_alpha_0", "alpha", [1, 2]), _region("t_beta_0", "beta", [3, 4])]
        row = design_row("t", d, regions, per_pair=True)
        _approx(float(row["mean_interdomain_pae"]), 10.0)  # 0.5*(12+8)
        _approx(float(row["max_interdomain_pae"]), 10.0)
        assert row["n_domains"] == 2
        _approx(float(row["pae_t_alpha_0_t_beta_0"]), 10.0)


def test_extract_dir_end_to_end_with_detections_json():
    import json

    with tempfile.TemporaryDirectory() as d:
        structs = os.path.join(d, "structs")
        pae_dir = os.path.join(d, "structs_pae")
        os.makedirs(structs)
        os.makedirs(pae_dir)

        # Two designs: "multi" (2 domains, scored), "single" (1 domain, N/A).
        for stem in ("multi", "single"):
            with open(os.path.join(structs, stem + ".pdb"), "w") as fh:
                fh.write("REMARK placeholder\nEND\n")

        multi_pae = np.full((4, 4), 2.0)
        multi_pae[0:2, 2:4] = 15.0
        multi_pae[2:4, 0:2] = 11.0
        _write_npz(os.path.join(pae_dir, "multi_pae.npz"), multi_pae, [1, 2, 3, 4])
        _write_npz(os.path.join(pae_dir, "single_pae.npz"), np.ones((3, 3)), [1, 2, 3])

        # EE detection JSON sidecar (so no EnzymeExplorer needed).
        detections = {
            "multi": [
                _region("multi_alpha_0", "alpha", [1, 2]),
                _region("multi_beta_0", "beta", [3, 4]),
            ],
            "single": [_region("single_alpha_0", "alpha", [1, 2, 3])],
        }
        json_path = os.path.join(d, "detections.json")
        with open(json_path, "w") as fh:
            json.dump(detections, fh)

        df = extract_interdomain_pae_dir(structs, pae_dir, detections_json=json_path)

        assert list(df.columns)[: len(BASE_COLUMNS)] == BASE_COLUMNS
        assert set(df["ID"]) == {"multi", "single"}
        assert os.path.isfile(structs + "_interdomain_pae.csv")

        m = df.set_index("ID").loc["multi"]
        _approx(float(m["mean_interdomain_pae"]), 13.0)  # 0.5*(15+11)
        assert int(m["n_domains"]) == 2

        s = df.set_index("ID").loc["single"]
        assert pd.isna(s["mean_interdomain_pae"])
        assert int(s["n_domains"]) == 1


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
