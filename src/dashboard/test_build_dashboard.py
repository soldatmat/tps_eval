from __future__ import annotations

"""Self-contained tests for dashboard/build_dashboard.py data-prep helpers.

Run from this directory (flat-module import, like the script does):
    cd src/dashboard && python test_build_dashboard.py
or under pytest:
    cd src/dashboard && python -m pytest test_build_dashboard.py -q

Exercises ONLY the pure/deterministic data-assembly logic (no HTML rendering,
no template): filename->tool label mapping, --designs spec parsing, CSV path
expansion, number sniffing, JSON sanitization, band-column compaction, metric
ordering, design-only grouping, the design-batch loader (synthetic temp CSVs),
and the large-mode base64 typed-array encoder. All inputs are synthetic and
in-memory / in a temp dir.
"""

import base64
import os
import struct
import tempfile

import build_dashboard as B


def test_parse_tool_label():
    assert B._parse_tool_label("/x/gen_plddt.csv") == "plddt"
    # longest-suffix-wins: the _self sibling keeps its own label
    assert B._parse_tool_label("/x/g_max_sequence_identity_self.csv") == "max_sequence_identity_self"
    # explicit remaps
    assert B._parse_tool_label("/x/g_motifs.csv") == "motif_search"
    assert (
        B._parse_tool_label("/x/g_embedding_esm1b_min_embedding_distance.csv")
        == "min_embedding_distance"
    )
    # unknown suffix -> full stem
    assert B._parse_tool_label("/x/weird_thing.csv") == "weird_thing"


def test_parse_design_specs():
    specs = B.parse_design_specs(["mine=a.csv,b.csv", "c.csv", "  ", "n2= d.csv "])
    assert specs == [
        ("mine", ["a.csv", "b.csv"]),
        (None, ["c.csv"]),
        ("n2", ["d.csv"]),
    ]


def test_parse_design_specs_equals_only_before_first_comma():
    # '=' after the first comma must NOT be treated as a name delimiter.
    specs = B.parse_design_specs(["a.csv,b=c.csv"])
    assert specs == [(None, ["a.csv", "b=c.csv"])]


def test_is_number():
    assert B._is_number("3.2") and B._is_number("5") and B._is_number(-1)
    assert not B._is_number("x") and not B._is_number(None) and not B._is_number("")


def test_ordered_metrics_known_first_extras_sorted():
    got = B._ordered_metrics(["zzz", "plddt", "aggregation", "aaa"])
    assert got == ["plddt", "aggregation", "aaa", "zzz"]


def test_sanitize_for_json_replaces_nonfinite():
    out = B._sanitize_for_json({"a": float("nan"), "b": [1.0, float("inf")], "c": "s", "d": 3})
    assert out == {"a": None, "b": [1.0, None], "c": "s", "d": 3}


def test_compact_numeric_and_categorical():
    col = {"mean": 1, "median": 2, "junk": 9, "std": 0.5}
    assert B._compact_numeric(col) == {"mean": 1, "median": 2, "std": 0.5}
    cat = {"count": 10, "n_missing": 1, "n_unique": 3, "frequencies": {"a": 5}, "junk": 1}
    out = B._compact_categorical(cat)
    assert out == {"count": 10, "n_missing": 1, "n_unique": 3, "frequencies": {"a": 5}}


def test_compact_column_carries_kind_and_by_strata():
    col = {
        "kind": "numeric", "mean": 1.0, "median": 2.0,
        "by_substrate": {"mono": {"mean": 3.0, "median": 4.0}, "bad": 7},
        "by_first_cyclization": {},
    }
    out = B._compact_column(col)
    assert out["kind"] == "numeric"
    assert out["mean"] == 1.0
    assert "by" in out and "substrate" in out["by"]
    assert out["by"]["substrate"]["mono"] == {"mean": 3.0, "median": 4.0}
    # non-dict strata dropped; empty labeling not carried
    assert "first_cyclization" not in out["by"]


def test_resolve_csv_paths_dir_glob_file_dedup():
    with tempfile.TemporaryDirectory() as d:
        for fn in ("a_plddt.csv", "b_soluprot.csv", "notes.txt"):
            open(os.path.join(d, fn), "w").close()
        # directory expands to the two csvs (sorted), .txt ignored
        got = B._resolve_csv_paths([d])
        assert [os.path.basename(p) for p in got] == ["a_plddt.csv", "b_soluprot.csv"]
        # explicit file + same dir -> deduped
        got2 = B._resolve_csv_paths([os.path.join(d, "a_plddt.csv"), d])
        assert sum(p.endswith("a_plddt.csv") for p in got2) == 1


def test_design_only_groups_skips_banded_and_orders():
    union_cols = {
        "mean_plddt": ("plddt", "numeric"),      # banded -> skipped
        "soluble": ("soluprot", "numeric"),      # not banded -> grouped
        "weird": (None, "categorical"),          # no tool -> "design metrics"
    }
    band_kind = {"mean_plddt": "numeric"}
    groups = B._design_only_groups(union_cols, band_kind)
    labels = [m for m, _ in groups]
    assert "plddt" not in labels
    assert "soluprot" in labels and "design metrics" in labels
    sol = dict(groups)["soluprot"]
    assert sol["soluble"] == {"kind": "numeric", "band_missing": True}


def test_load_design_batch_merges_and_types():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "gen_plddt.csv")
        with open(p, "w") as f:
            f.write("ID,mean_plddt,domain_architecture\n")
            f.write("x,90.5,alpha-beta\n")
            f.write("y,,gamma\n")       # missing numeric -> None
            f.write("z,70,alpha\n")
        ds = B.load_design_batch([p], {}, name="t")
        assert ds is not None
        assert ds["n"] == 3 and ds["ids"] == ["x", "y", "z"]
        assert ds["col_kind"]["mean_plddt"] == "numeric"
        assert ds["col_kind"]["domain_architecture"] == "categorical"
        assert ds["values"]["mean_plddt"] == [90.5, None, 70.0]
        assert ds["col_tool"]["mean_plddt"] == "plddt"


def test_load_design_batch_self_columns_suffixed():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "gen_local_sequence_search_self.csv")
        with open(p, "w") as f:
            f.write("ID,local_sequence_identity\nx,42\n")
        ds = B.load_design_batch([p], {}, name="t")
        # self-file columns get a _self suffix so they don't clobber gen-vs-train.
        assert "local_sequence_identity_self" in ds["values"]
        assert "local_sequence_identity" not in ds["values"]


def test_load_design_batch_none_when_no_csv():
    with tempfile.TemporaryDirectory() as d:
        assert B.load_design_batch([d], {}, name="t") is None


def test_load_design_batch_skips_wide_feature_matrix():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "gen_embedding.csv")
        ncols = B._MAX_DESIGN_COLUMNS_PER_FILE + 5
        with open(p, "w") as f:
            f.write("ID," + ",".join(f"d{i}" for i in range(ncols)) + "\n")
            f.write("x," + ",".join("0.1" for _ in range(ncols)) + "\n")
        # only that (too-wide) file present -> nothing usable -> None
        assert B.load_design_batch([p], {}, name="t") is None


def test_encode_large_values_roundtrip():
    ds = {
        "col_kind": {"mp": "numeric", "arch": "categorical"},
        "values": {
            "mp": [90.5, None, 70.0],
            "arch": ["alpha", "beta", None],
        },
    }
    B._encode_large_values([ds])
    assert "values" not in ds
    enc = ds["values_enc"]
    # numeric -> float32, null -> NaN
    mp = enc["mp"]
    decoded = struct.unpack("<%df" % mp["n"], base64.b64decode(mp["b64"]))
    assert decoded[0] == 90.5 and decoded[2] == 70.0
    assert decoded[1] != decoded[1]  # NaN
    # categorical -> category table + uint16 codes, null -> 0xFFFF
    arch = enc["arch"]
    codes = struct.unpack("<%dH" % arch["n"], base64.b64decode(arch["b64"]))
    assert arch["cats"] == ["alpha", "beta"]
    assert codes == (0, 1, B._NULL_CODE)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
