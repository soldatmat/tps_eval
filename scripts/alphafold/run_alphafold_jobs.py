import os
import subprocess
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Run AlphaFold jobs from CSV")
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file')
    parser.add_argument('--working_directory', required=True, help='Working directory for all runs')
    parser.add_argument('--save_directory', required=False, default=None, help='Directory to save the structures (default: working_directory/structs)')
    parser.add_argument('--id_column_name', required=False, default='ID', help='Name of the column containing sequence IDs to be used as names (default: ID)')
    parser.add_argument('--sequence_column_name', required=False, default='sequence', help='Name of the column containing sequences (default: sequence)')
    parser.add_argument('--cluster', required=False, default='aurum', help='Cluster on which this script is run. (default: aurum)')
    parser.add_argument('--submit_args', required=False, default=None, help='Additional arguments to pass to the cluster job submission')
    args = parser.parse_args()

    df = pd.read_csv(args.csv_path)
    working_directory = str(args.working_directory)

    run_alphafold_jobs(
        df,
        working_directory,
        args.id_column_name,
        args.sequence_column_name,
        args.cluster,
        save_directory=args.save_directory,
        submit_args=args.submit_args,
    )


def run_alphafold_jobs(
    df,
    working_directory,
    id_column_name,
    sequence_column_name,
    cluster,
    save_directory=None,
    submit_args=None,
):
    n_skipped = 0
    job_ids = []
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
        source_path = os.path.realpath(__file__)
        script_dir = os.path.dirname(source_path)
        cmd = [
            'bash',
            f"{script_dir}/submit_job.sh",
            '--cluster', cluster,
            '--job_name', "alphafold",
            '--job_args', f"{working_directory} {sequence_id} {sequence} {save_directory if save_directory else ''}",
        ]
        if submit_args:
            cmd += ['--submit_args', submit_args]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        cmd_output = result.stdout.strip()
        job_id = cmd_output.split()[-1] # Last word printed by `submit_job.sh` is the job ID
        job_ids.append(job_id)
    print(f"Skipped running AlphaFold for {n_skipped} sequences without Uniprot ID that already have PDB files.")
    
    print("AlphaFold job IDs:", job_ids)
    # Last thing printed has to be list of the job IDs for job dependencies.

if __name__ == "__main__":
    main()
