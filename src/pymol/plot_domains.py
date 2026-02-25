import argparse
from pathlib import Path
import pandas as pd
import pickle
from pymol import cmd
from tqdm import tqdm

from src.pymol.utils import show_organic_and_metals
from src.pymol.constants import PNG_WIDTH, PNG_DPI, PNG_RAY, SMALL_SECONDARY_STRUCTURE_WIDTH


def parse_args() -> argparse.Namespace:
    """
    This function parses arguments
    :return: current argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="A script to compare sequence similarity of new protein structures to selected known structures.")
    parser.add_argument("--structures_selection_csv", type=str, required=True, help="Path to CSV file with structure selection.")
    parser.add_argument("--structures_column_name", type=str, default="ID", help="Name of the column in the CSV file that contains structure paths.")
    parser.add_argument("--structures_file_suffix", type=str, default="", help="Suffix to add to structure names when looking for structure files (e.g. if structure files are named like 'structure_suffix.pdb').")
    parser.add_argument("--structures_root", type=str, required=True, help="Path to new structures")
    parser.add_argument("--domain_structures_root", type=str, required=True, help="Path to new structures")
    parser.add_argument("--domains_pkl", type=str, required=True, help="Path to domains pickle file")
    parser.add_argument("--output_root", type=str, required=True, help="Path to output the images and PyMOL sessions.")
    parser.add_argument("--store_pymol_sessions", action=argparse.BooleanOptionalAction, default=True, help="Flag to store PyMOL sessions for each comparison.")
    parser.add_argument("--rerun_existing", action=argparse.BooleanOptionalAction, default=True, help="Flag to rerun the script even if output already exists.")
    return parser.parse_args()


def load_enzyme_explorer_domains(domains_pkl_path: str) -> dict:
    class DummyClass:
        def __init__(self, *args, **kwargs):
            pass

    class UnknownClassUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            try:
                # Try to import normally
                return super().find_class(module, name)
            except (ModuleNotFoundError, AttributeError):
                # Replace unknown classes with dict
                return DummyClass

    def to_dict(obj):
        if isinstance(obj, dict):
            return {k: to_dict(v) for k, v in obj.items()}
        elif hasattr(obj, "__dict__"):
            return {k: to_dict(v) for k, v in obj.__dict__.items()}
        elif isinstance(obj, (list, tuple, set)):
            return type(obj)(to_dict(x) for x in obj)
        else:
            return obj

    with open(domains_pkl_path, "rb") as f:
        data = UnknownClassUnpickler(f).load()
    clean_data = to_dict(data)
    return clean_data


def main(args: argparse.Namespace):
    n_skipped_existing = 0
    n_skipped_missing_structures = 0
    skipped_missing_structures = []

    with open(args.structures_selection_csv) as f:
        df = pd.read_csv(f)
        filename_2_known_regions = load_enzyme_explorer_domains(args.domains_pkl)

        progress_bar = tqdm(df.iterrows(), total=len(df), bar_format="{l_bar}{bar}| {percentage:3.0f}%")
        for _, row in progress_bar:
            # Extract structure names and paths
            tqdm.write("")
            structure_name = row[args.structures_column_name]
            structure_file_name = structure_name + args.structures_file_suffix
            structure_path = Path(args.structures_root) / f"{structure_file_name}.pdb"
            progress_bar.set_description(f"Processing {structure_file_name:<{50}} ")
            tqdm.write(f"Processing '{structure_file_name}' ...")

            # Prepare output paths and check if output already exists
            output_dir = Path(args.output_root) / structure_file_name
            session_output_path = output_dir / "domains.pse"
            domains_output_path = output_dir / "domains.png"
            if (not args.rerun_existing) and domains_output_path.exists():
                tqdm.write(f"Output for structure '{structure_file_name}' already exists at '{domains_output_path.parent}'. Skipping '{structure_file_name}'...")
                n_skipped_existing += 1
                continue

            # Check if known regions are available and if structure file exists
            if not structure_path.exists():
                tqdm.write(f"Warning: Structure file '{structure_path}' does not exist. Skipping '{structure_file_name}'...")
                n_skipped_missing_structures += 1
                skipped_missing_structures.append(structure_file_name)
                continue
            if structure_name not in filename_2_known_regions:
                tqdm.write(f"Warning: No known regions found for structure '{structure_file_name}' in domains pickle file. Skipping '{structure_file_name}'...")
                n_skipped_missing_structures += 1
                skipped_missing_structures.append(structure_file_name)
                continue
            domain_names = [domain['module_id'] for domain in filename_2_known_regions[structure_name]]
            if len(domain_names) != len(set(domain_names)):
                tqdm.write(f"Warning: Duplicate domain names found for structure '{structure_file_name}'. Skipping '{structure_file_name}'...")
                n_skipped_missing_structures += 1
                skipped_missing_structures.append(structure_file_name)
                continue
            missing_domain_file = False
            domain_structure_paths = []
            for domain_name in domain_names:
                domain_structure_path = Path(args.domain_structures_root) / f"{domain_name}.pdb"
                if not domain_structure_path.exists():
                    tqdm.write(f"Warning: Domain structure file '{domain_structure_path}' does not exist. Skipping '{structure_file_name}'...")
                    n_skipped_missing_structures += 1
                    skipped_missing_structures.append(structure_file_name)
                    missing_domain_file = True
                    break
                domain_structure_paths.append(domain_structure_path)
            if missing_domain_file:
                continue

            # PyMOL: structure with detected domains
            cmd.reinitialize()
            cmd.load(str(structure_path), "structure")
            for i, domain_name in enumerate(domain_names):
                cmd.load(str(domain_structure_paths[i]), domain_name)
            cmd.set("cartoon_oval_width", SMALL_SECONDARY_STRUCTURE_WIDTH, "structure")
            cmd.set("cartoon_rect_width", SMALL_SECONDARY_STRUCTURE_WIDTH, "structure")
            show_organic_and_metals(cmd)
            cmd.orient()
            output_dir.mkdir(parents=True, exist_ok=True)
            if args.store_pymol_sessions:
                cmd.save(str(session_output_path))
            cmd.png(str(domains_output_path), width=PNG_WIDTH, dpi=PNG_DPI, ray=PNG_RAY)
            tqdm.write(f"Saved image for '{structure_file_name}' at '{domains_output_path.parent}'")

    print(f"\nFinished processing. Skipped {n_skipped_existing} runs with existing output and {n_skipped_missing_structures} runs with missing structures.")
    if skipped_missing_structures:
        print("Skipped runs with missing structures:")
        for structure_file_name in skipped_missing_structures:
            print(f" • {structure_file_name}")


if __name__ == "__main__":
    args = parse_args()
    main(args)
