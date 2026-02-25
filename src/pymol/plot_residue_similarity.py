"""
Takes a list of new structures and a list of preselected matches among known structures
and colors the new structures by their sequence similarity to the known ones.

The similarity is defined as the BLOSUM90 score of aligned residues, calculated by the
color_by_mutation function from the color_by_mutations.py script.
"""

import argparse
from pathlib import Path
import pandas as pd
from pymol import cmd
from tqdm import tqdm

from vendor.pymol_scripts.color_by_mutations import color_by_mutation


def parse_args() -> argparse.Namespace:
    """
    This function parses arguments
    :return: current argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="A script to compare sequence similarity of new protein structures to selected known structures.")
    parser.add_argument("--structures_selection_csv", type=str, required=True, help="Path to CSV file with structure selection.")
    parser.add_argument("--structures_column_name", type=str, default="query", help="Name of the column in the CSV file that contains structure paths.")
    parser.add_argument("--known_structures_column_name", type=str, default="max_alntmscore_target", help="Name of the column in the CSV file that contains paths to known structures to compare to.")
    parser.add_argument("--structures_root", type=str, required=True, help="Path to new structures")
    parser.add_argument("--known_structures_root", type=str, required=True, help="Directory containing structures of known proteins")
    parser.add_argument("--output_root", type=str, required=True, help="Path to output the images and PyMOL sessions.")
    parser.add_argument("--store_pymol_sessions", action=argparse.BooleanOptionalAction, default=True, help="Flag to store PyMOL sessions for each comparison.")
    parser.add_argument("--rerun_existing", action=argparse.BooleanOptionalAction, default=True, help="Flag to rerun the script even if output already exists.")
    return parser.parse_args()


def main(args: argparse.Namespace):
    n_skipped_existing = 0
    n_skipped_missing_structures = 0
    skipped_missing_structures = []

    with open(args.structures_selection_csv) as f:
        df = pd.read_csv(f)
        progress_bar = tqdm(df.iterrows(), total=len(df), bar_format="{l_bar}{bar}| {percentage:3.0f}%")
        for _, row in progress_bar:
            # Extract structure names and paths
            tqdm.write("")
            structure_name = row[args.structures_column_name]
            if structure_name.endswith("_A"):
                tqdm.write(f"Warning: Removing '_A' suffix from structure name '{structure_name}'")
                structure_name = structure_name[:-2]
            structure_path = Path(args.structures_root) / f"{structure_name}.pdb"
            known_structure_name = row[args.known_structures_column_name]
            known_structure_path = Path(args.known_structures_root) / f"{known_structure_name}.pdb"
            run_name = f"{structure_name}-{known_structure_name}"
            progress_bar.set_description(f"Processing {run_name:<{60}} ")
            tqdm.write(f"Processing {run_name} ...")

            # Prepare output paths and check if output already exists
            output_dir = Path(args.output_root) / run_name
            session_output_path = output_dir / "color_by_mutation.pse"
            color_by_mutation_output_path = output_dir / "color_by_mutation.png"
            alignment_output_path = output_dir / "alignment.png"
            if (not args.rerun_existing) and color_by_mutation_output_path.exists() and alignment_output_path.exists():
                tqdm.write(f"Output for structure '{structure_name}' already exists at '{color_by_mutation_output_path.parent}'. Skipping '{run_name}'...")
                n_skipped_existing += 1
                continue

            # Check if structure files exist
            if not structure_path.exists():
                tqdm.write(f"Error: Structure file '{structure_path}' does not exist. Skipping '{run_name}'...")
                n_skipped_missing_structures += 1
                skipped_missing_structures.append(run_name)
                continue
            if not known_structure_path.exists():
                tqdm.write(f"Error: Known structure file '{known_structure_path}' does not exist. Skipping '{run_name}'...")
                n_skipped_missing_structures += 1
                skipped_missing_structures.append(run_name)
                continue

            # PyMOL: color_by_mutation
            cmd.reinitialize()
            cmd.load(str(structure_path), "structure")
            cmd.load(str(known_structure_path), "known")
            color_by_mutation("structure", "known", verbosity=0)
            cmd.hide("everything", "known")
            cmd.show("sticks", "structure and organic")
            cmd.color("atomic", "structure and organic")
            cmd.color("gray", "structure and organic and elem C")
            cmd.show("spheres", "structure and metals")
            cmd.orient()
            output_dir.mkdir(parents=True, exist_ok=True)
            if args.store_pymol_sessions:
                cmd.save(str(session_output_path))
            cmd.png(str(color_by_mutation_output_path), width=2000, dpi=300, ray=1)

            # PyMOL: alignment
            cmd.show("cartoon", "known and polymer")
            cmd.color("sulfur", "structure and polymer")
            cmd.color("skyblue", "known and polymer")
            cmd.png(str(alignment_output_path), width=2000, dpi=300, ray=1)
            tqdm.write(f"Saved images for '{run_name}' at '{color_by_mutation_output_path.parent}'")

    print(f"\nFinished processing. Skipped {n_skipped_existing} runs with existing output and {n_skipped_missing_structures} runs with missing structures.")
    if skipped_missing_structures:
        print("Skipped runs with missing structures:")
        for run_name in skipped_missing_structures:
            print(f" • {run_name}")

if __name__ == "__main__":
    args = parse_args()
    main(args)
