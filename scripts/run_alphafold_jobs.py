import os
import subprocess
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Run AlphaFold jobs from CSV")
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file')
    parser.add_argument('--working_directory', required=True, help='Working directory for all runs')
    parser.add_argument('--id_column_name', required=False, default='ID', help='Name of the column containing sequence IDs to be used as names (default: ID)')
    parser.add_argument('--sequence_column_name', required=False, default='sequence', help='Name of the column containing sequences (default: sequence)')
    parser.add_argument('--cluster', required=False, default='aurum', help='Cluster on which this script is run. (default: aurum)')
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    working_directory = str(args.working_directory)

    run_alphafold_jobs(df, working_directory, args.id_column_name, args.sequence_column_name, args.cluster)


def run_alphafold_jobs(df, working_directory, id_column_name, sequence_column_name, cluster):
    n_skipped = 0
    for _, row in df.iterrows():
        # if pd.notna(row.get('Uniprot_ID')):
        #     continue

        sequence_id = row[id_column_name]

        output_pdb_path = os.path.join(working_directory, "structs", f"{sequence_id}.pdb")
        if os.path.exists(output_pdb_path):
            n_skipped += 1
            continue
        
        sequence = row[sequence_column_name]

        print(f"Running AlphaFold for sequence ID: {sequence_id}")
        cmd = [
            'bash',
            f"./{cluster}/run_alphafold.sh",
            '--working_directory', working_directory,
            '--sequence_id', str(sequence_id),
            '--sequence', str(sequence)
        ]
        subprocess.run(cmd, check=True)
    print(f"Skipped running AlphaFold for {n_skipped} sequences without Uniprot ID that already have PDB files.")


if __name__ == "__main__":
    main()
