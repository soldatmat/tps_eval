import pandas as pd
import argparse
import subprocess

import sys
import os
import importlib.util

def main():
    """
    Run AlphaFold predictions for sequences listed in a CSV file.

    The CSV must contain a column named "Uniprot_ID" (the downloader expects such a column).
    This column is used to check for precomputed AlphaFold structures.
    Generated input & output files are named using the column specified by --id_column_name.
    """
    parser = argparse.ArgumentParser(description="Run AlphaFold from CSV")
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file')
    parser.add_argument('--working_directory', required=True, help='Working directory for all runs')
    parser.add_argument('--id_column_name', required=False, default='ID', help='Name of the column containing sequence IDs to be used as names (default: ID)')
    parser.add_argument('--sequence_column_name', required=False, default='sequence', help='Name of the column containing sequences (default: sequence)')
    parser.add_argument('--cluster', required=False, default='aurum', help='Cluster on which this script is run. (default: aurum)')
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    working_directory = str(args.working_directory)

    # Run alphafold_struct_downloader to download structures for sequences with UniProt IDs
    spec = importlib.util.spec_from_file_location("alphafold_struct_downloader", "./alphafold_struct_downloader.py")
    alphafold_struct_downloader = importlib.util.module_from_spec(spec)
    sys.modules["alphafold_struct_downloader"] = alphafold_struct_downloader
    spec.loader.exec_module(alphafold_struct_downloader)
    structs_output_path = working_directory + "/structs"
    alphafold_struct_downloader.main(
        structures_output_path=structs_output_path,
        path_to_file_with_ids=args.csv_path,
        n_jobs=1,
        save_name_column_name=args.id_column_name,
    )

    # Run a separate AlphaFold job for each sequence in the CSV without UniProt ID
    n_skipped = 0
    for _, row in df.iterrows():
        # if pd.notna(row.get('Uniprot_ID')):
        #     continue

        sequence_id = row[args.id_column_name]
        sequence = row[args.sequence_column_name]

        output_pdb_path = os.path.join(working_directory, "structs", f"{sequence_id}.pdb")
        if os.path.exists(output_pdb_path):
            n_skipped += 1
            continue

        print(f"Running AlphaFold for sequence ID: {sequence_id}")
        cmd = [
            'bash',
            f"./{args.cluster}/run_alphafold.sh",
            '--working_directory', working_directory,
            '--sequence_id', str(sequence_id),
            '--sequence', str(sequence)
        ]
        subprocess.run(cmd, check=True)
    print(f"Skipped running AlphaFold for {n_skipped} sequences without Uniprot ID that already have PDB files.")

if __name__ == "__main__":
    main()
