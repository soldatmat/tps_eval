import argparse

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_cif", required=True, type=str)
    parser.add_argument("--output_pdb", required=True, type=str)
    args = parser.parse_args()

    cif_parser = MMCIFParser()
    structure = cif_parser.get_structure("structure", args.input_cif)

    io = PDBIO()
    io.set_structure(structure)
    io.save(args.output_pdb, write_end=True)


if __name__ == "__main__":
    main()
