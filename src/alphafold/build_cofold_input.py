"""Build the AlphaFold3 fan-out input CSV(s) for a given co-fold mode.

Used by ``scripts/run_alphafold_fanout.sh``. Turns a generated FASTA into one or more
``run_alphafold_jobs.py`` input CSVs plus a manifest, encoding the active-site co-fold:

  none            apo (protein only)
  mg              + trinuclear Mg2+ cluster (3x CCD ``MG``)
  mg_ppi          + Mg cluster + bare diphosphate head group (CCD ``POP``)
  mg_<substrate>  + Mg cluster + ONE forced prenyl-PP substrate (SMILES) for EVERY design
                  (<substrate> in gpp|fpp|ggpp|gfpp)
  mg_ee           + Mg cluster + each design's OWN EnzymeExplorer-predicted substrate
                  (requires --enzymeexplorer_csv; designs whose EE argmax is not co-foldable fall back
                  to Mg-only)

Ions are CCD-coded (passed to run_alphafold_jobs as ``--ion_*``); substrates are SMILES
(passed as ``--ligand_*``). ``mg_ee`` cannot share one ligand column (substrate varies per
design and run_alphafold_jobs does not skip empty cells), so it is split into one CSV per
substrate group + a Mg-only group; the manifest tells the caller which CSVs carry a ligand.

Output: writes ``<out_dir>/af3_input[_<tag>].csv`` files and ``<out_dir>/af3_cofold_manifest.tsv``
(columns: csv_path, has_ligand, n_designs). Column names are fixed: ions ``ion{i}_id``/
``ion{i}_ccd``; ligand ``lig1_id``/``lig1_smiles``.
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from alphafold.cofold_substrates import SUBSTRATE_SMILES, COFOLDABLE, is_cofoldable  # noqa: E402

ION_ID_COLS = ["ion1_id", "ion2_id", "ion3_id", "ion4_id"]
ION_CCD_COLS = ["ion1_ccd", "ion2_ccd", "ion3_ccd", "ion4_ccd"]
LIG_ID_COL = "lig1_id"
LIG_SMILES_COL = "lig1_smiles"

VALID_MODES = ["none", "mg", "mg_ppi", "mg_ee"] + [f"mg_{s.lower()}" for s in COFOLDABLE]


def ions_for(cofold: str) -> List[Tuple[str, str]]:
    """(ion_id, CCD) list. 3x Mg for every holo mode; + POP only for mg_ppi (the substrate
    modes carry their own diphosphate via the ligand)."""
    if cofold == "none":
        return []
    ions = [("MG1", "MG"), ("MG2", "MG"), ("MG3", "MG")]
    if cofold == "mg_ppi":
        ions.append(("PPI", "POP"))
    return ions


def read_fasta(path: str) -> List[Tuple[str, str]]:
    recs, cur, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    recs.append((cur, "".join(seq)))
                cur, seq = line[1:].split()[0], []
            elif line.strip():
                seq.append(line.strip())
    if cur is not None:
        recs.append((cur, "".join(seq)))
    return recs


def _write_csv(path: str, recs: List[Tuple[str, str]], ions: List[Tuple[str, str]],
               ligand_smiles: Optional[str]) -> None:
    id_cols = ION_ID_COLS[:len(ions)]
    ccd_cols = ION_CCD_COLS[:len(ions)]
    header = ["ID", "sequence"]
    for idc, cdc in zip(id_cols, ccd_cols):
        header += [idc, cdc]
    if ligand_smiles is not None:
        header += [LIG_ID_COL, LIG_SMILES_COL]
    ion_vals = [v for (iid, ccd) in ions for v in (iid, ccd)]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for rid, s in recs:
            row = [rid, s] + ion_vals
            if ligand_smiles is not None:
                row += ["LIG", ligand_smiles]
            w.writerow(row)


def _ee_substrate_per_design(ee_csv: str) -> Dict[str, str]:
    """{ID: substrate_code} from the EE seq-only argmax (reusing the substrate_class loader)."""
    from knn.substrate_class import load_ee_substrate  # noqa: E402
    return {cid: code for cid, (code, _score) in load_ee_substrate(ee_csv).items()}


def build(fasta: str, cofold: str, out_dir: str, ee_csv: Optional[str] = None) -> List[Tuple[str, bool, int]]:
    """Write the job CSV(s); return list of (csv_path, has_ligand, n_designs)."""
    if cofold not in VALID_MODES:
        raise SystemExit(f"Unknown --cofold {cofold!r}; valid: {', '.join(VALID_MODES)}")
    os.makedirs(out_dir, exist_ok=True)
    recs = read_fasta(fasta)
    ions = ions_for(cofold)
    manifest: List[Tuple[str, bool, int]] = []

    if cofold in ("none", "mg", "mg_ppi"):
        p = os.path.join(out_dir, "af3_input.csv")
        _write_csv(p, recs, ions, ligand_smiles=None)
        manifest.append((p, False, len(recs)))

    elif cofold.startswith("mg_") and cofold[3:].upper() in SUBSTRATE_SMILES:
        code = cofold[3:].upper()
        p = os.path.join(out_dir, "af3_input.csv")
        _write_csv(p, recs, ions, ligand_smiles=SUBSTRATE_SMILES[code])
        manifest.append((p, True, len(recs)))

    elif cofold == "mg_ee":
        if not ee_csv:
            raise SystemExit("--cofold mg_ee requires --enzymeexplorer_csv (the EnzymeExplorer seq-only "
                             "CSV); run enzyme_explorer_sequence_only first.")
        codes = _ee_substrate_per_design(ee_csv)
        groups: Dict[str, List[Tuple[str, str]]] = {}
        for rid, s in recs:
            code = codes.get(rid, "")
            key = code.upper() if is_cofoldable(code) else "_MGONLY"
            groups.setdefault(key, []).append((rid, s))
        for key in sorted(groups):
            grp = groups[key]
            if key == "_MGONLY":
                p = os.path.join(out_dir, "af3_input_mgonly.csv")
                _write_csv(p, grp, ions, ligand_smiles=None)
                manifest.append((p, False, len(grp)))
            else:
                p = os.path.join(out_dir, f"af3_input_{key.lower()}.csv")
                _write_csv(p, grp, ions, ligand_smiles=SUBSTRATE_SMILES[key])
                manifest.append((p, True, len(grp)))

    man_path = os.path.join(out_dir, "af3_cofold_manifest.tsv")
    with open(man_path, "w") as fh:
        fh.write("csv_path\thas_ligand\tn_designs\n")
        for p, has_lig, n in manifest:
            fh.write(f"{p}\t{int(has_lig)}\t{n}\n")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description="Build AF3 fan-out input CSV(s) for a co-fold mode.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--cofold", required=True, choices=VALID_MODES)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--enzymeexplorer_csv", default=None, help="EE seq-only CSV (required for mg_ee).")
    args = ap.parse_args()
    manifest = build(args.fasta, args.cofold, args.output_dir, ee_csv=args.enzymeexplorer_csv)
    n_ion = len(ions_for(args.cofold))
    # Quoted so the caller can `eval` these lines safely (values are space-separated lists).
    print(f'ION_ID_COLS="{" ".join(ION_ID_COLS[:n_ion])}"')
    print(f'ION_CCD_COLS="{" ".join(ION_CCD_COLS[:n_ion])}"')
    print(f'LIG_ID_COL="{LIG_ID_COL}"')
    print(f'LIG_SMILES_COL="{LIG_SMILES_COL}"')
    for p, has_lig, n in manifest:
        print(f"[build_cofold_input] {os.path.basename(p)}: {n} design(s), "
              f"ligand={'yes' if has_lig else 'no'}")


if __name__ == "__main__":
    main()
