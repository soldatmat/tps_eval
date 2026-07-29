"""Unit tests for domain_composition.py (TPS domain composition table builder).

Run: python test_domain_composition.py   (no pytest / EnzymeExplorer needed).

Tests only the pure, EE-free logic: regions_to_rows (per-type counts, ordered
architecture string, zero-domain fill, unexpected-type handling), the JSON sidecar
parser, and the default save-path derivation. Detection itself (detect_domains_json)
imports EnzymeExplorer lazily and is NOT exercised here.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile


from tps_eval.enzyme_explorer.domain_composition import (
    COLUMNS,
    _default_save_path,
    load_detections_json,
    regions_to_rows,
)


def test_regions_to_rows_counts_and_architecture():
    seq_to_regions = {
        "s1": [
            {"domain": "beta", "module_id": "s1_beta_0"},
            {"domain": "alpha", "module_id": "s1_alpha_1"},
            {"domain": "alpha", "module_id": "s1_alpha_0"},
        ],
    }
    df = regions_to_rows(seq_to_regions, ["s1", "s2"]).set_index("ID")
    assert list(df.reset_index().columns) == COLUMNS
    # s1: 3 domains, 2 alpha, 1 beta; architecture sorted by (type, index).
    assert int(df.loc["s1", "n_domains"]) == 3
    assert int(df.loc["s1", "n_alpha"]) == 2
    assert int(df.loc["s1", "n_beta"]) == 1
    assert df.loc["s1", "domain_architecture"] == "alpha-alpha-beta"
    # s2: absent from detections -> zero-domain row.
    assert int(df.loc["s2", "n_domains"]) == 0
    assert df.loc["s2", "domain_architecture"] == ""
    print("ok regions_to_rows_counts_and_architecture")


def test_unexpected_domain_type_counted_without_column():
    seq_to_regions = {"s1": [{"domain": "weird", "module_id": "s1_weird_0"},
                             {"domain": "alpha", "module_id": "s1_alpha_0"}]}
    df = regions_to_rows(seq_to_regions, ["s1"]).set_index("ID")
    # Unexpected type still counts toward n_domains + architecture but has no n_* column.
    assert int(df.loc["s1", "n_domains"]) == 2
    assert "weird" in df.loc["s1", "domain_architecture"]
    assert "n_weird" not in df.columns
    print("ok unexpected_domain_type_counted_without_column")


def test_load_detections_json_roundtrip():
    tmp = tempfile.mkdtemp(prefix="domcomp_")
    p = os.path.join(tmp, "det.json")
    payload = {"s1": [{"module_id": "s1_alpha_0", "domain": "alpha", "tmscore": 0.9}]}
    with open(p, "w") as fh:
        json.dump(payload, fh)
    parsed = load_detections_json(p)
    assert parsed["s1"][0]["domain"] == "alpha"
    print("ok load_detections_json_roundtrip")


def test_default_save_path():
    assert _default_save_path("/a/b/structs").endswith("structs_domain_composition.csv")
    # Trailing separator handled.
    assert _default_save_path("/a/b/structs/").endswith("structs_domain_composition.csv")
    print("ok default_save_path")


def _install_fake_enzymeexplorer(calls, *, raises: bool = False):
    """Stub out EE so detect_domains_json runs without EnzymeExplorer installed.

    The fake ``detect_domains`` mimics the two things that made the real one
    collide across concurrent tools: it creates each per-domain subdirectory with
    a bare ``mkdir(parents=True)`` (no ``exist_ok`` — EE's own check-then-mkdir),
    and it writes the secondary-structure pickle. Returns a restore callable.
    """
    import types

    saved = {k: v for k, v in sys.modules.items() if k.split(".")[0] in ("enzymeexplorer", "yaml")}

    def _fake_detect_domains(args):
        calls.append(
            {
                "domains_output_path": args.domains_output_path,
                "secondary_structure_residues_path": args.secondary_structure_residues_path,
                "store_domains": args.store_domains,
            }
        )
        # EE creates these unconditionally-ish; a shared path across two runs blows up.
        for domain in ("alpha", "gamma"):
            os.makedirs(os.path.join(args.domains_output_path, domain))
        with open(args.secondary_structure_residues_path, "wb") as fh:
            fh.write(b"ssr")
        if raises:
            raise RuntimeError("boom")
        with open(str(args.detections_output_path)[: -len(".pkl")] + ".json", "w") as fh:
            json.dump({"s1": [{"module_id": "s1_alpha_0", "domain": "alpha"}]}, fh)

    pkg = types.ModuleType("enzymeexplorer")
    src = types.ModuleType("enzymeexplorer.src")
    sp = types.ModuleType("enzymeexplorer.src.structure_processing")
    dd = types.ModuleType("enzymeexplorer.src.structure_processing.domain_detections")
    dd.DEFAULT_DOMAIN_TEMPLATES = [{"name": "alpha"}]
    dd.detect_domains = _fake_detect_domains
    for name, mod in [
        ("enzymeexplorer", pkg),
        ("enzymeexplorer.src", src),
        ("enzymeexplorer.src.structure_processing", sp),
        ("enzymeexplorer.src.structure_processing.domain_detections", dd),
    ]:
        sys.modules[name] = mod
    try:  # PyYAML is only guaranteed inside the EE env
        import yaml  # noqa: F401
    except ImportError:
        fake_yaml = types.ModuleType("yaml")
        fake_yaml.safe_dump = lambda obj: str(obj)
        sys.modules["yaml"] = fake_yaml

    def restore():
        for name in list(sys.modules):
            if name.split(".")[0] in ("enzymeexplorer", "yaml"):
                del sys.modules[name]
        sys.modules.update(saved)

    return restore


def test_detect_domains_json_scratch_is_per_invocation_and_cleaned():
    """Two runs sharing an output dir must not share EE scratch paths.

    Regression: domain_composition and interdomain_pae run concurrently on the
    same structures dir, and both used a fixed `_ee_domains_scratch` /
    `_ee_secondary_structure_residues.pkl` sibling of their JSON — the loser died
    with FileExistsError on `_ee_domains_scratch/gamma`.
    """
    from tps_eval.enzyme_explorer.domain_composition import detect_domains_json

    tmp = tempfile.mkdtemp(prefix="domcomp_scratch_")
    calls = []
    restore = _install_fake_enzymeexplorer(calls)
    try:
        for stem in ("structs_domain_composition_detections", "structs_interdomain_pae_detections"):
            out = detect_domains_json(tmp, out_json_path=os.path.join(tmp, stem + ".json"))
            assert out["s1"][0]["domain"] == "alpha"
    finally:
        restore()

    assert len(calls) == 2
    scratch_dirs = [c["domains_output_path"] for c in calls]
    ssr_paths = [c["secondary_structure_residues_path"] for c in calls]
    assert scratch_dirs[0] != scratch_dirs[1], scratch_dirs
    assert ssr_paths[0] != ssr_paths[1], ssr_paths
    # Nothing consumes the per-domain PDBs, so EE is told not to write them.
    assert all(c["store_domains"] is False for c in calls)
    # Scratch is removed once detection returns — no leftovers in the output dir.
    for path in scratch_dirs + ssr_paths:
        assert not os.path.exists(path), path
    leftovers = [n for n in os.listdir(tmp) if "_ee_scratch_" in n]
    assert not leftovers, leftovers
    print("ok detect_domains_json_scratch_is_per_invocation_and_cleaned")


def test_detect_domains_json_cleans_scratch_on_failure():
    from tps_eval.enzyme_explorer.domain_composition import detect_domains_json

    tmp = tempfile.mkdtemp(prefix="domcomp_scratch_fail_")
    calls = []
    restore = _install_fake_enzymeexplorer(calls, raises=True)
    try:
        try:
            detect_domains_json(tmp, out_json_path=os.path.join(tmp, "d.json"))
            raise AssertionError("expected the stubbed detect_domains to raise")
        except RuntimeError:
            pass
    finally:
        restore()

    assert not os.path.exists(calls[0]["domains_output_path"])
    assert not [n for n in os.listdir(tmp) if "_ee_scratch_" in n]
    print("ok detect_domains_json_cleans_scratch_on_failure")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
