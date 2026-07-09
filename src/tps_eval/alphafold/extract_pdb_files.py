"""Convert AlphaFold3 mmCIF model(s) to PDB, preserving the B-factor (pLDDT) and the
ligand/ion HETATMs, with PDB-INVALID residue names sanitized.

Why this exists (not vendor/cif_to_pdb): AF3 holo folds carry the co-folded ligand chain
through the mmCIF with a comp_id like ``LIG_B`` (>3 chars). PDB's resName is a fixed 3-char
field, so writing ``LIG_B`` overflows into the chainID column and the resulting .pdb is
unparseable — Biopython's PDBParser dies with ``invalid literal for int() with base 10: 'B'``,
which made every downstream structure tool (plddt, active_site_geometry, ion_site_check,
substrate_positioning, pocket_descriptors, ...) read 0 residues and emit empty rows. We
truncate any >3-char resName to a valid 3-char code (``LIG_B`` -> ``LIG``) before writing.
Ions / standard residues (``MG``, ``POP``, ``MET`` — all <=3 chars) are untouched, so
ion_site_check and the composition-based substrate detection in substrate_positioning still
work. (The B-factor pLDDT preservation that motivated vendor/cif_to_pdb is identical here:
Biopython MMCIFParser reads B_iso into atom.bfactor and PDBIO writes it back.)
"""
import argparse
import re
import warnings
from pathlib import Path

from Bio.PDB import PDBIO, Select
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class _ResidueSelect(Select):
    """protein/nucleic -> ATOM; ligand/ion/cofactor -> HETATM (kept); water -> dropped.
    (Biopython hetero flag is residue.id[0]: ' ' standard, 'W' water, 'H_*' hetero.)"""

    def __init__(self, keep_hetero: bool = True):
        self.keep_hetero = keep_hetero

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        if hetflag == " ":
            return True
        if hetflag == "W":
            return False
        return self.keep_hetero


def _sanitize_resnames(structure) -> None:
    """Truncate any residue resName longer than the PDB 3-char field to a valid 3-char
    alphanumeric code (e.g. AF3's 'LIG_B' -> 'LIG'). In place."""
    for residue in structure.get_residues():
        name = residue.resname.strip()
        if len(name) > 3:
            residue.resname = (re.sub(r"[^A-Za-z0-9]", "", name)[:3].upper() or "LIG")


def cif_to_pdb_sanitized(input_cif: str, output_pdb: str, keep_hetero: bool = True) -> None:
    """mmCIF -> PDB (Biopython): keep the full complex (protein ATOM + ligand/ion HETATM),
    preserve the B-factor pLDDT, drop water, and sanitize PDB-invalid resNames so the output
    parses. Pass keep_hetero=False for a protein-only structure."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = MMCIFParser(QUIET=True).get_structure("s", input_cif)
    _sanitize_resnames(structure)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, _ResidueSelect(keep_hetero=keep_hetero))


def main():
    parser = argparse.ArgumentParser(
        description="Convert AF3 mmCIF -> PDB (B-factor/pLDDT preserved, ligand HETATMs kept, "
        "long resNames sanitized). Single-file mode (--input_cif/--output_pdb) or directory "
        "mode (--af_output_folder/--struct_folder).")
    parser.add_argument("--input_cif", type=str, help="Single mmCIF to convert.")
    parser.add_argument("--output_pdb", type=str, help="Output PDB path (with --input_cif).")
    parser.add_argument("--af_output_folder", type=str,
                        help="AF3 af_output dir with per-structure subfolders.")
    parser.add_argument("--struct_folder", type=str,
                        help="Output dir for extracted PDBs (with --af_output_folder).")
    parser.add_argument("--no_hetero", dest="keep_hetero", action="store_false",
                        help="Protein/nucleic atoms only; drop ligands/ions.")
    parser.set_defaults(keep_hetero=True)
    args = parser.parse_args()

    if args.input_cif and args.output_pdb:
        cif_to_pdb_sanitized(args.input_cif, args.output_pdb, keep_hetero=args.keep_hetero)
    elif args.af_output_folder and args.struct_folder:
        Path(args.struct_folder).mkdir(parents=True, exist_ok=True)
        for folder in Path(args.af_output_folder).iterdir():
            input_cif = folder / (folder.name + "_model.cif")
            if not input_cif.is_file():
                continue
            output_pdb = Path(args.struct_folder) / (folder.name + ".pdb")
            cif_to_pdb_sanitized(str(input_cif), str(output_pdb), keep_hetero=args.keep_hetero)
    else:
        parser.error("give either --input_cif + --output_pdb, or --af_output_folder + --struct_folder")


if __name__ == "__main__":
    main()
