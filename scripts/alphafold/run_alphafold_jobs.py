import os
import subprocess
import argparse
import pandas as pd


ID_SEQUENCE_SEPARATOR = ' '
SEQUENCE_PAIR_SEPARATOR = ' '
SEED_SEPARATOR = ' '


def parse_args():
    parser = argparse.ArgumentParser(description="Run AlphaFold jobs from CSV")

    # CSV file
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file')
    parser.add_argument('--protein_id_column_names', type=str, nargs='+', default=['ID'], help='List of protein id column names.')
    parser.add_argument('--protein_sequence_column_names', type=str, nargs='+', default=['sequence'], help='List of protein sequence column names in the same order as --protein_id_column_names.')
    parser.add_argument('--ligand_id_column_names', required=False, type=str, nargs='+', help='List containing ligand id column names.')
    parser.add_argument('--ligand_smiles_column_names', required=False, type=str, nargs='+', help='List containing ligand SMILES column names in the same order as --ligand_id_column_names.')
    parser.add_argument('--csv_delimiter', type=str, default=',', help='Delimiter used in the CSV file (default: ,)')

    # AlphaFold parameters
    parser.add_argument('--model_seeds', type=int, nargs='+', default=[42], help='List of model seeds to use for AlphaFold runs.')

    # Save directory
    parser.add_argument('--working_directory', required=True, help='Working directory for all runs')
    parser.add_argument('--save_directory', default=None, help='Directory to save the structures (default: working_directory/structs)')
    
    # Job submission scripts
    parser.add_argument('--cluster', default='aurum', help='Cluster on which this script is run. (default: aurum)')
    parser.add_argument('--submit_args', type=str, default="", help='Additional arguments to pass to the cluster job submission')
    parser.add_argument("--skip_existing", default=True, action=argparse.BooleanOptionalAction, help='Skip sequences that already have PDB files generated.')
    
    args = parser.parse_args()
    return args


def prepare_submit_args(submit_args, cluster, default_job_name, working_directory):
    if cluster == 'aurum':
        if "--job-name=" not in submit_args:
            submit_args += f" --job-name=AF_{default_job_name}"
        if "--output=" not in submit_args:
            submit_args += f" --output={working_directory}/logs/%x.%j.out"

    submit_args = f'"{submit_args.strip()}"'

    return submit_args


def run_alphafold_jobs(
    df,
    working_directory,
    cluster,
    protein_column_names,
    ligand_column_names=[],
    save_directory=None,
    submit_args="",
    model_seeds=[42],
    skip_existing=True,
):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    model_seeds = SEED_SEPARATOR.join([str(seed) for seed in model_seeds])

    n_skipped = 0
    job_ids = []
    for _, row in df.iterrows():
        proteins = [
            (row[protein_id_column_name], row[protein_sequence_column_name])
            for protein_id_column_name, protein_sequence_column_name in protein_column_names
        ]
        ligands = [
            (row[ligand_id_column_name], row[ligand_smiles_column_name])
            for ligand_id_column_name, ligand_smiles_column_name in ligand_column_names
        ]
        combined_protein_ids = "_".join([id for id, sequence in proteins])
        proteins = SEQUENCE_PAIR_SEPARATOR.join([f'{id}{ID_SEQUENCE_SEPARATOR}{sequence}' for id, sequence in proteins])
        ligands = SEQUENCE_PAIR_SEPARATOR.join([f'{id}{ID_SEQUENCE_SEPARATOR}{smiles}' for id, smiles in ligands])

        output_pdb_path = os.path.join(working_directory, "structs", f"{combined_protein_ids}.pdb")
        if skip_existing and os.path.exists(output_pdb_path):
            n_skipped += 1
            continue

        print(f"Running AlphaFold for protein ID: {combined_protein_ids}")
        job_args = f"\"--working_directory {working_directory} --sequence_id {combined_protein_ids} --proteins {proteins} --ligands {ligands}{f" --save_directory {save_directory}" if save_directory else ''} --model_seeds {model_seeds}\""
        job_submit_args = prepare_submit_args(submit_args, cluster=cluster, default_job_name=combined_protein_ids, working_directory=working_directory)
        cmd = [
            'bash',
            f"{script_dir}/../submit_job.sh",
            '--cluster', cluster,
            '--job_name', "alphafold",
            '--job_args', job_args,
            '--submit_args', job_submit_args,
        ]

        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error submitting AlphaFold job for protein ID: {combined_protein_ids}")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            print("exit code:", result.returncode)
            exit(1)
        # print(result.stdout) # debug

        cmd_output = result.stdout.strip()
        job_id = cmd_output.split()[-1] # Last word printed by `submit_job.sh` is the job ID
        job_ids.append(job_id)
    print(f"Skipped running AlphaFold for {n_skipped} sequences without Uniprot ID that already have PDB files.")
    
    print("AlphaFold job IDs:", job_ids)
    # Last thing printed has to be list of the job IDs for job dependencies.

def main():
    args = parse_args()
    df = pd.read_csv(args.csv_path, sep=args.csv_delimiter)

    protein_column_names = list(zip(args.protein_id_column_names, args.protein_sequence_column_names))
    ligand_column_names = []
    if args.ligand_id_column_names and args.ligand_smiles_column_names:
        ligand_column_names = list(zip(args.ligand_id_column_names, args.ligand_smiles_column_names))

    run_alphafold_jobs(
        df,
        args.working_directory,
        args.cluster,
        protein_column_names,
        ligand_column_names,
        save_directory=args.save_directory,
        submit_args=args.submit_args,
        model_seeds=args.model_seeds,
        skip_existing=args.skip_existing,
    )


if __name__ == "__main__":
    main()
