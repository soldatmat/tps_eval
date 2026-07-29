from __future__ import annotations

"""Self-contained tests for cofold_substrates.py (the class-I TPS prenyl-diphosphate
SMILES lookup table used for AF3 holo co-folding).

Run from this directory (so the package-style imports inside the module resolve):
    cd src/alphafold && python test_cofold_substrates.py
or under pytest:
    cd src/alphafold && python -m pytest test_cofold_substrates.py -q

Pure hardcoded data — no model, no I/O. Locks in table integrity (COFOLDABLE ==
the SMILES keys, uppercase vocabulary, no accidental duplicates via case, every
SMILES carries a diphosphate head + balanced parens) and the case-insensitive
accessor / membership contract that build_cofold_input relies on.
"""

import os
import sys
from pathlib import Path


from tps_eval.alphafold.cofold_substrates import (
    COFOLDABLE,
    SUBSTRATE_SMILES,
    is_cofoldable,
    smiles_for,
)


def test_cofoldable_is_a_subset_of_the_table_keys():
    """Every co-foldable code must have a SMILES — but NOT the reverse.

    These two were once identical, and this test used to assert equality. They are
    deliberately different now: ``SUBSTRATE_SMILES`` is the full substrate vocabulary
    shared with CataPro, while ``COFOLDABLE`` is the narrower AF3 co-folding policy.
    DMAPP (C5) and C35 have SMILES for CataPro but are not co-folded — see
    ``is_cofoldable``'s docstring, which names exactly those two.
    """
    missing = set(COFOLDABLE) - set(SUBSTRATE_SMILES)
    assert not missing, f"co-foldable codes with no SMILES: {missing}"
    # No dup entries in the ordered list.
    assert len(COFOLDABLE) == len(set(COFOLDABLE))
    # Lock in the intended asymmetry so a silent re-merge is caught.
    assert set(SUBSTRATE_SMILES) - set(COFOLDABLE) == {"DMAPP", "C35"}
    for code in ("DMAPP", "C35"):
        assert not is_cofoldable(code)
        assert smiles_for(code)  # has a SMILES for CataPro


def test_codes_uppercase_and_unique_case_insensitive():
    """Vocabulary codes are UPPERCASE and unique even folded to lower case."""
    for code in SUBSTRATE_SMILES:
        assert code == code.upper(), code
    lowered = [c.lower() for c in SUBSTRATE_SMILES]
    assert len(lowered) == len(set(lowered)), "case-insensitive duplicate code"


def test_every_smiles_has_diphosphate_and_balanced_parens():
    """Each entry is a prenyl-diphosphate: carries the OP(...)OP(...) head and has
    balanced parentheses (a cheap SMILES sanity gate, no RDKit)."""
    for code, smi in SUBSTRATE_SMILES.items():
        assert smi and isinstance(smi, str), code
        assert "OP(" in smi, f"{code}: no diphosphate head in {smi!r}"
        assert smi.count("(") == smi.count(")"), f"{code}: unbalanced parens"
        assert smi.count("[") == smi.count("]"), f"{code}: unbalanced brackets"


def test_chain_length_monotonic():
    """GPP<FPP<GGPP<GFPP: SMILES length grows with the prenyl chain (C10<C15<C20<C25)."""
    lengths = [len(SUBSTRATE_SMILES[c]) for c in ["GPP", "FPP", "GGPP", "GFPP"]]
    assert lengths == sorted(lengths), lengths
    assert len(set(lengths)) == 4, "distinct chain lengths expected"


def test_smiles_for_case_insensitive():
    assert smiles_for("gpp") == SUBSTRATE_SMILES["GPP"]
    assert smiles_for("GPP") == SUBSTRATE_SMILES["GPP"]
    assert smiles_for("Fpp") == SUBSTRATE_SMILES["FPP"]


def test_smiles_for_unknown_raises():
    for bad in ("edsq", "cpp", "2xfpp", "xyz"):
        try:
            smiles_for(bad)
        except KeyError:
            continue
        raise AssertionError(f"expected KeyError for non-cofoldable code {bad!r}")


def test_is_cofoldable():
    assert is_cofoldable("gpp")
    assert is_cofoldable("GGPP")
    assert not is_cofoldable("")
    assert not is_cofoldable("EDSQ")
    assert not is_cofoldable("CPP")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"All {len(tests)} tests passed.")


if __name__ == "__main__":
    main()
