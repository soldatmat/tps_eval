"""
Adapted from original script from Raman Samusevich

This script downloads protein structures predicted by AlphaFold2
"""

import argparse
import logging
from pathlib import Path
from functools import partial

import requests
from multiprocessing import Pool

import pandas as pd

from tqdm.auto import tqdm

logging.basicConfig()
logger = logging.getLogger("Downloading AlphaFold2 structures")
logger.setLevel(logging.INFO)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def parse_args() -> argparse.Namespace:
    """
    This function parses arguments
    :return: current argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--structures_output_path", type=str, default="../data/af_structs"
    )
    parser.add_argument("--path_to_file_with_ids", type=str, 
                        default="../data/uniprot_ids_of_interest.txt", help="Path to a file containing UniProt IDs,"
                                                                            "for which the script will download AF2 structures")
    parser.add_argument("--n_jobs", type=int, default=1)

    args = parser.parse_args()
    return args

def download_af_struct(uniprot_id, root_af, fails_count=0, max_fails_count=3):
    if isinstance(uniprot_id, tuple):
        uniprot_id, save_name = uniprot_id
    else:
        save_name = uniprot_id

    try:
        URL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb"
        response = requests.get(URL)
        with open(root_af / f"{save_name}.pdb", "wb") as file:
            file.write(response.content)
    except:
        logger.warning(f"Error downloading AlphaFold2 structure for {uniprot_id}")
        if fails_count < max_fails_count:
            download_af_struct(uniprot_id, root_af, fails_count+1)


def main(
        structures_output_path,
        path_to_file_with_ids,
        n_jobs,
        uniprot_id_column_name='Uniprot_ID',
        save_name_column_name='Enzyme_marts_ID',
    ):
    """
    This function downloads protein structures predicted by AlphaFold
    """
    root_af = Path(structures_output_path)
    if not root_af.exists():
        root_af.mkdir()

    download_af_struct_for_current_root = partial(download_af_struct, root_af=root_af)

    if path_to_file_with_ids.endswith('.txt'):
        with open(path_to_file_with_ids, 'r') as file:
            all_ids_of_interest = [line.strip() for line in file.readlines()]
    elif path_to_file_with_ids.endswith('.csv'):
        df = pd.read_csv(path_to_file_with_ids)
        if uniprot_id_column_name not in df.columns:
            print(f"Column with Uniprot IDs '{uniprot_id_column_name}' not found in the CSV file. No structures downloaded.")
            return
        df = df.dropna(subset=[uniprot_id_column_name])
        all_ids_of_interest = df[uniprot_id_column_name].astype(str).tolist()
        save_names = df[save_name_column_name].astype(str).tolist()
        all_ids_of_interest = list(zip(all_ids_of_interest, save_names))
    else:
        raise ValueError("Unsupported file format for UniProt IDs. Please provide a .txt or .csv file.")
    
    if all_ids_of_interest == []:
        logger.info("No UniProt IDs found in the provided file. No structures downloaded.")
        return
    
    # Filter out IDs for which the structure file already exists
    if isinstance(all_ids_of_interest[0], tuple):
        filtered_ids = [
            (uniprot_id, save_name)
            for uniprot_id, save_name in all_ids_of_interest
            if not (root_af / f"{save_name}.pdb").exists()
        ]
    else:
        filtered_ids = [
            uniprot_id
            for uniprot_id in all_ids_of_interest
            if not (root_af / f"{uniprot_id}.pdb").exists()
        ]
    skipped_count = len(all_ids_of_interest) - len(filtered_ids)
    logger.info(f"Skipped downloading structures for {skipped_count} sequences because the PDB file already exists.")

    with Pool(processes=n_jobs) as pool:
        pool.map(download_af_struct_for_current_root, filtered_ids)


if __name__ == "__main__":
    args = parse_args()
    main(args.structures_output_path, args.path_to_file_with_ids, args.n_jobs)
