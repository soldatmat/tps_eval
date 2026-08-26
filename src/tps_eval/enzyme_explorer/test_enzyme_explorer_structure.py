"""Tests for the structure-mode EnzymeExplorer reshaper.

The behaviour under test that actually matters for the gated metric: a design
with no detected TPS domain must come out as a FAILURE (NaN structure score),
never silently backfilled with the PLM-only fallback score.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tps_eval.enzyme_explorer.enzyme_explorer_structure import (
    DOMAINS_CSV,
    FALLBACK_CSV,
    STATUS_MISSING,
    STATUS_NO_DOMAINS,
    STATUS_SCORED,
    build_structure_predictions,
    default_save_path,
    design_ids_from_fasta,
    load_id_map,
    reshape_to_csv,
    summarize,
)

# EnzymeExplorer's `main`-branch schema.
MAIN_COLS = ["id", "TPS_p", "GPP_p", "FPP_p", "GGPP_p", "GFPP_p", "CPP_p", "EDSQ_p",
             "2xFPP_p", "2xGGPP_p", "IDS_p", "TPS_raw", "GPP_raw", "FPP_raw", "GGPP_raw",
             "GFPP_raw", "CPP_raw", "EDSQ_raw", "2xFPP_raw", "2xGGPP_raw", "IDS_raw",
             "sequence"]
# The older `revision`-branch schema.
REVISION_COLS = ["id", "sequence", "TPS_score", "GPP_score", "TPS_p_calibrated",
                 "GPP_p_calibrated"]


def _write(path, columns, rows):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _main_row(design_id, tps_raw, tps_cal):
    row = {c: 0.1 for c in MAIN_COLS}
    row.update(id=design_id, sequence="MAAA", TPS_raw=tps_raw, TPS_p=tps_cal)
    return row


def test_no_domains_design_fails_and_is_not_backfilled(tmp_path):
    """The core contract: no domain hit -> NaN structure score, fallback kept apart."""
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.91, 0.95)])
    _write(tmp_path / FALLBACK_CSV, MAIN_COLS, [_main_row("d2", 0.88, 0.93)])

    df = build_structure_predictions(str(tmp_path), ["d1", "d2"]).set_index("ID")

    assert df.loc["d1", "ee_struct_status"] == STATUS_SCORED
    assert df.loc["d1", "ee_struct_domain_hit"] is True or bool(df.loc["d1", "ee_struct_domain_hit"])
    assert df.loc["d1", "isTPS_struct_raw"] == pytest.approx(0.91)

    # d2 scored 0.88 by the FALLBACK classifier -- a high score that must NOT
    # become a structure score, or the gate is fooled exactly as before.
    assert df.loc["d2", "ee_struct_status"] == STATUS_NO_DOMAINS
    assert not bool(df.loc["d2", "ee_struct_domain_hit"])
    assert pd.isna(df.loc["d2", "isTPS_struct_raw"])
    assert pd.isna(df.loc["d2", "isTPS_struct_cal"])
    assert df.loc["d2", "isTPS_fallback_raw"] == pytest.approx(0.88)

    # ...and the NaN fails a >= test, which is how the gate consumes it.
    assert not (df["isTPS_struct_raw"] >= 0.5)["d2"]
    assert (df["isTPS_struct_raw"] >= 0.5)["d1"]


def test_design_absent_from_both_csvs_is_missing_not_dropped(tmp_path):
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.91, 0.95)])
    _write(tmp_path / FALLBACK_CSV, MAIN_COLS, [_main_row("d2", 0.10, 0.10)])

    df = build_structure_predictions(str(tmp_path), ["d1", "d2", "d3_never_folded"])

    assert len(df) == 3
    row = df.set_index("ID").loc["d3_never_folded"]
    assert row["ee_struct_status"] == STATUS_MISSING
    assert pd.isna(row["isTPS_struct_raw"])
    assert pd.isna(row["isTPS_fallback_raw"])


def test_header_only_domains_csv_means_every_design_rejected(tmp_path):
    """The observed real case: 56/56 composition artifacts, zero domain hits.

    EE writes an empty domains frame as a header-only file carrying just
    ``id,sequence`` -- it must read as "nobody scored", not as a parse error.
    """
    _write(tmp_path / DOMAINS_CSV, ["id", "sequence"], [])
    _write(tmp_path / FALLBACK_CSV, MAIN_COLS,
           [_main_row(f"d{i}", 0.7, 0.8) for i in range(5)])

    df = build_structure_predictions(str(tmp_path), [f"d{i}" for i in range(5)])

    assert (df["ee_struct_status"] == STATUS_NO_DOMAINS).all()
    assert df["isTPS_struct_raw"].isna().all()
    assert (df["isTPS_struct_raw"] >= 0.5).sum() == 0
    assert summarize(df) == {
        "n_designs": 5, "n_domain_hit": 0, "n_no_domains": 5,
        "n_missing": 0, "n_isTPS_struct_pass": 0,
    }


def test_revision_branch_schema_is_accepted(tmp_path):
    """`<class>_score` / `<class>_p_calibrated` resolve like `_raw` / `_p`."""
    _write(tmp_path / DOMAINS_CSV, REVISION_COLS,
           [{"id": "d1", "sequence": "MAAA", "TPS_score": 0.77, "GPP_score": 0.2,
             "TPS_p_calibrated": 0.85, "GPP_p_calibrated": 0.3}])
    _write(tmp_path / FALLBACK_CSV, ["id", "sequence"], [])

    df = build_structure_predictions(str(tmp_path), ["d1"]).set_index("ID")

    assert df.loc["d1", "isTPS_struct_raw"] == pytest.approx(0.77)
    assert df.loc["d1", "isTPS_struct_cal"] == pytest.approx(0.85)
    assert df.loc["d1", "GPP_struct_raw"] == pytest.approx(0.2)


def test_overlapping_ids_raise(tmp_path):
    """The two passes partition the input; an overlap means mismatched CSVs."""
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.9, 0.9)])
    _write(tmp_path / FALLBACK_CSV, MAIN_COLS, [_main_row("d1", 0.4, 0.4)])

    with pytest.raises(ValueError, match="BOTH"):
        build_structure_predictions(str(tmp_path), ["d1"])


def test_duplicate_design_ids_raise(tmp_path):
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.9, 0.9)])
    _write(tmp_path / FALLBACK_CSV, ["id", "sequence"], [])

    with pytest.raises(ValueError, match="duplicate design IDs"):
        build_structure_predictions(str(tmp_path), ["d1", "d1"])


def test_missing_output_dir_yields_all_missing(tmp_path):
    """A prediction that never ran must fail every design, not crash."""
    df = build_structure_predictions(str(tmp_path / "absent"), ["d1", "d2"])
    assert (df["ee_struct_status"] == STATUS_MISSING).all()
    assert (df["isTPS_struct_raw"] >= 0.5).sum() == 0


def test_design_universe_from_fasta_and_default_path(tmp_path):
    fasta = tmp_path / "GEN.fasta"
    fasta.write_text(">d1 some description\nMAAA\n>d2\nMBBB\n")
    assert design_ids_from_fasta(str(fasta)) == ["d1", "d2"]

    out_dir = tmp_path / "GEN_enzyme_explorer"
    out_dir.mkdir()
    _write(out_dir / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.91, 0.95)])
    _write(out_dir / FALLBACK_CSV, MAIN_COLS, [_main_row("d2", 0.2, 0.2)])

    structs = tmp_path / "GEN_esmfold_structs"
    structs.mkdir()
    assert default_save_path(str(structs)).endswith(
        "GEN_esmfold_structs_enzyme_explorer_structure.csv"
    )

    path = reshape_to_csv(str(out_dir), structs_dir=str(structs), fasta_path=str(fasta))
    df = pd.read_csv(path)
    assert list(df["ID"]) == ["d1", "d2"]


def test_structs_dir_fallback_universe_misses_unfolded_designs(tmp_path):
    """Documented limitation of the .pdb-stem universe, pinned by a test."""
    structs = tmp_path / "GEN_esmfold_structs"
    structs.mkdir()
    (structs / "d1.pdb").write_text("ATOM\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _write(out_dir / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.91, 0.95)])
    _write(out_dir / FALLBACK_CSV, ["id", "sequence"], [])

    path = reshape_to_csv(str(out_dir), structs_dir=str(structs))
    df = pd.read_csv(path)
    assert list(df["ID"]) == ["d1"]  # d2, which never folded, is simply not seen


def test_id_map_restores_original_ids(tmp_path):
    """A sanitized run must come back keyed by the ORIGINAL design IDs."""
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write(out_dir / DOMAINS_CSV, MAIN_COLS, [_main_row("cell_s0_5_00001", 0.9, 0.95)])
    _write(out_dir / FALLBACK_CSV, ["id", "sequence"], [])
    id_map = tmp_path / "id_map.csv"
    pd.DataFrame({"sanitized_id": ["cell_s0_5_00001"],
                  "original_id": ["cell_s0.5_00001"]}).to_csv(id_map, index=False)
    structs = tmp_path / "GEN_esmfold_structs"; structs.mkdir()
    fasta = tmp_path / "GEN.fasta"; fasta.write_text(">cell_s0_5_00001\nMAAA\n")

    path = reshape_to_csv(str(out_dir), structs_dir=str(structs),
                         fasta_path=str(fasta), id_map_csv=str(id_map))
    assert list(pd.read_csv(path)["ID"]) == ["cell_s0.5_00001"]


def test_incomplete_id_map_refuses(tmp_path):
    """A partial map would key the table by two ID conventions -- refuse instead."""
    out_dir = tmp_path / "out"; out_dir.mkdir()
    _write(out_dir / DOMAINS_CSV, MAIN_COLS, [_main_row("a_1", 0.9, 0.95)])
    _write(out_dir / FALLBACK_CSV, MAIN_COLS, [_main_row("b_2", 0.1, 0.1)])
    id_map = tmp_path / "id_map.csv"
    pd.DataFrame({"sanitized_id": ["a_1"], "original_id": ["a.1"]}).to_csv(id_map, index=False)
    structs = tmp_path / "S"; structs.mkdir()
    fasta = tmp_path / "GEN.fasta"; fasta.write_text(">a_1\nMAAA\n>b_2\nMBBB\n")

    with pytest.raises(ValueError, match="does not cover"):
        reshape_to_csv(str(out_dir), structs_dir=str(structs), fasta_path=str(fasta),
                       id_map_csv=str(id_map))


def test_id_map_collision_refused(tmp_path):
    id_map = tmp_path / "id_map.csv"
    pd.DataFrame({"sanitized_id": ["a_1", "a_2"], "original_id": ["a.1", "a.1"]}
                 ).to_csv(id_map, index=False)
    with pytest.raises(ValueError, match="collided"):
        load_id_map(str(id_map))


def test_one_byte_csv_is_an_empty_branch_not_an_error(tmp_path):
    """EE writes a bare newline when a frame has no columns at all.

    This is the COMMON case (every design got a domain hit -> empty fallback),
    and pd.read_csv raises EmptyDataError on it. Reading it as "nobody fell back"
    is the only correct behaviour; raising would abort an otherwise-fine arm.
    """
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS,
           [_main_row(f"d{i}", 0.95, 0.99) for i in range(3)])
    (tmp_path / FALLBACK_CSV).write_text("\n")  # exactly what EE writes

    df = build_structure_predictions(str(tmp_path), [f"d{i}" for i in range(3)])

    assert (df["ee_struct_status"] == STATUS_SCORED).all()
    assert df["ee_struct_domain_hit"].all()
    assert df["isTPS_fallback_raw"].isna().all()
    assert (df["isTPS_struct_raw"] >= 0.5).sum() == 3


def test_zero_byte_csv_is_also_tolerated(tmp_path):
    _write(tmp_path / DOMAINS_CSV, MAIN_COLS, [_main_row("d1", 0.95, 0.99)])
    (tmp_path / FALLBACK_CSV).write_bytes(b"")

    df = build_structure_predictions(str(tmp_path), ["d1"])
    assert df.loc[0, "ee_struct_status"] == STATUS_SCORED


def test_both_csvs_empty_means_every_design_missing(tmp_path):
    """A prediction that produced nothing at all must fail every design, not crash."""
    (tmp_path / DOMAINS_CSV).write_text("\n")
    (tmp_path / FALLBACK_CSV).write_text("\n")

    df = build_structure_predictions(str(tmp_path), ["d1", "d2"])
    assert (df["ee_struct_status"] == STATUS_MISSING).all()
    assert (df["isTPS_struct_raw"] >= 0.5).sum() == 0
