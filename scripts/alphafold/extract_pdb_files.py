import argparse
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO

from src.alphafold.cif_to_pdb import main as cif_to_pdb_main


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--af_output_folder", required=True, type=str, help="AlphaFold output folder containing subfolders for each structure.")
    parser.add_argument("--struct_folder", required=True, type=str, help="Folder where the extraceted PDBs will be saved.")
    args = parser.parse_args()

    for folder in Path(args.af_output_folder).iterdir():
        input_cif = folder / (folder.name + "_model.cif")
        output_pdb = Path(args.struct_folder) / (folder.name + ".pdb")

        cif_to_pdb_main_args = argparse.Namespace(input_cif=str(input_cif), output_pdb=str(output_pdb))
        cif_to_pdb_main(cif_to_pdb_main_args)


if __name__ == "__main__":
    main()
