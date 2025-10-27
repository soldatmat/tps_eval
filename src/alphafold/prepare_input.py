import argparse
import json


SEQUENCE_IDS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence_id", type=str, default=None)
    parser.add_argument("--sequence", type=str, default=None)
    parser.add_argument("--proteins", type=str, nargs="+", default=None, help="List of proteins in format: ID1 SEQ1 ID2 SEQ2 ... All following tokens (until next --option) are parsed as proteins.")
    parser.add_argument("--ligands", type=str, nargs="+", default=[], help="List of ligands in format: ID1 SMILES1 ID2 SMILES2 ... All following tokens (until next --option) are parsed as ligands.")
    parser.add_argument("--save_path", required=True, type=str)
    parser.add_argument("--model_seeds", type=int, nargs="+", default=[42])
    args = parser.parse_args()

    if (args.sequence_id is None) ^ (args.sequence is None):
        parser.error("Both --sequence_id and --sequence must be provided together.")
    if not (args.proteins is not None or (args.sequence_id is not None and args.sequence is not None)):
        parser.error("Either provide --proteins or both --sequence_id and --sequence.")

    return args


def pair_ids_and_sequences(sequences_arg):
    sequences_pairs = []
    for i in range(0, len(sequences_arg), 2):
        id = sequences_arg[i]
        seq = sequences_arg[i + 1]
        sequences_pairs.append((id, seq))

    return sequences_pairs


def format_data(args):
    if len(args.proteins) + len(args.ligands) > len(SEQUENCE_IDS):
        raise ValueError(f"Number of proteins and ligands exceeds the number of predefined sequence ids ({len(args.proteins)} + {len(args.ligands)} > {len(SEQUENCE_IDS)}).")

    if args.sequence is not None:
        name = args.sequence_id
        proteins = [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": args.sequence,
                }
            },
        ]
    elif args.proteins is not None:
        protein_tuples = pair_ids_and_sequences(args.proteins)
        name = protein_tuples[0][0]
        proteins = [
            {
                "protein": {
                    "id": [SEQUENCE_IDS[index]],
                    "sequence": seq,
                }
            }
            for index, (id, seq) in enumerate(protein_tuples)
        ]
    else:
        raise ValueError("No proteins provided.")

    ligands = [
        {
            "ligand": {
                "id": [SEQUENCE_IDS[index + len(protein_tuples)]],
                "smiles": smiles,
            }
        }
        for index, (id, smiles) in enumerate(pair_ids_and_sequences(args.ligands))
    ]

    data = {
        "name": name,
        "modelSeeds": args.model_seeds,
        "sequences": [
                *proteins,
                *ligands,
        ],
        "dialect": "alphafold3",
        "version": 2,
    }

    return data


def main():
    args = parse_args()
    data = format_data(args)
    with open(args.save_path, "w") as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    main()
