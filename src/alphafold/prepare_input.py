import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence_id", type=str)
    parser.add_argument("--sequence", type=str)
    parser.add_argument("--save_path", type=str)
    args = parser.parse_args()

    data = {
        "name": args.sequence_id,
        "modelSeeds": [42],
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": args.sequence,
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 2,
    }

    with open(args.save_path, "w") as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    main()
