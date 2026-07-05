from __future__ import annotations

"""Self-contained tests for prepare_input.py (builds the AF3 input JSON dict).

Run from this directory:
    cd src/alphafold && python test_prepare_input.py
or under pytest:
    cd src/alphafold && python -m pytest test_prepare_input.py -q

Pure dict construction — no AF3, no I/O (we call format_data directly with an
argparse.Namespace). Locks in the chain-id assignment (proteins get A,B,C..., then
ligands, then ions off the shared SEQUENCE_IDS pool), the single-sequence path (which
regressed with a NoneType crash), and the ligand/ion sub-structure shape.
"""

import argparse
import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from alphafold.prepare_input import format_data, pair_ids_and_sequences


def _ns(**kw):
    base = dict(sequence_id=None, sequence=None, proteins=None, ligands=[], ions=[],
                model_seeds=[42], save_path="x")
    base.update(kw)
    return argparse.Namespace(**base)


def test_pair_ids_and_sequences():
    assert pair_ids_and_sequences(["a", "SEQA", "b", "SEQB"]) == [("a", "SEQA"), ("b", "SEQB")]
    assert pair_ids_and_sequences([]) == []


def test_single_sequence_path():
    """The --sequence / --sequence_id path (proteins=None) must not crash and must
    emit one protein on chain A. (Regression: format_data used len(args.proteins)
    unconditionally -> TypeError when proteins is None.)"""
    data = format_data(_ns(sequence_id="prot1", sequence="MKTAAR"))
    assert data["name"] == "prot1"
    assert data["dialect"] == "alphafold3" and data["version"] == 2
    assert data["modelSeeds"] == [42]
    assert data["sequences"] == [{"protein": {"id": ["A"], "sequence": "MKTAAR"}}]


def test_multi_protein_chain_ids():
    """Proteins get A, B, C ... in order; name defaults to first protein id."""
    data = format_data(_ns(proteins=["p1", "MK", "p2", "GG", "p3", "AA"]))
    assert data["name"] == "p1"
    ids = [entry["protein"]["id"][0] for entry in data["sequences"]]
    seqs = [entry["protein"]["sequence"] for entry in data["sequences"]]
    assert ids == ["A", "B", "C"]
    assert seqs == ["MK", "GG", "AA"]


def test_ligands_and_ions_get_following_chain_ids():
    """1 protein (A) + 2 ligands (B, C) + 1 ion (D): chain ids come off the shared
    SEQUENCE_IDS pool, offset by the counts of the earlier blocks."""
    data = format_data(_ns(
        proteins=["p1", "MK"],
        ligands=["l1", "CCO", "l2", "c1ccccc1"],
        ions=["i1", "MG"],
    ))
    seqs = data["sequences"]
    # order in output: proteins, ligands, ions
    assert seqs[0] == {"protein": {"id": ["A"], "sequence": "MK"}}
    assert seqs[1] == {"ligand": {"id": ["B"], "smiles": "CCO"}}
    assert seqs[2] == {"ligand": {"id": ["C"], "smiles": "c1ccccc1"}}
    assert seqs[3] == {"ligand": {"id": ["D"], "ccdCodes": ["MG"]}}


def test_sequence_id_overrides_name_for_proteins():
    data = format_data(_ns(sequence_id="custom", proteins=["p1", "MK"]))
    assert data["name"] == "custom"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
