# Pipeline for in silico evaluation of proteins on computational clusters
- Use `tps_eval/scripts/submit_job.sh` to submit jobs on any cluster.
- Every evaluation tool has its main run script created as `tps_eval/scripts/run_<tool_name>`.
- The main run script can be run from a job script on a computational cluster created as `tps_eval/scripts/<cluster>/<jobs>/<tool_name>.sh`.

# Installation
```sh
cd tps_eval
. setup.sh
```

# Running AlphaFold
Currently, Alphafold jobs are configured only for the IOCB Aurum cluster.

1) Install https://github.com/soldatmat/tps_eval

2) Run alphafold jobs like so:
```sh
cd "tps_eval/scripts/alphafold"
conda activate tps_eval
```
```sh
python run_alphafold_jobs.py \
    --csv_path /path/to/file_with_proteins_and_ligands.csv \
    --working_directory /path/to/directory/for/results/alphafold_structs
```
## Optional arguments
For usage, see `python run_alphafold_jobs.py -h`
### Example of each optional argument:
```sh
python run_alphafold_jobs.py \
    --csv_path /path/to/file_with_proteins_and_ligands.csv \
    --protein_id_column_names Protein_1 Protein_2 \
    --protein_sequence_column_names Sequence_1 Sequence_2 \
    --ligand_id_column_names Substrate Ligand \
    --ligand_smiles_column_names Substrate_SMILE Ligand_SMILE \
    --csv_delimiter ';' \
    --model_seeds 42 101 \
    --working_directory /path/to/directory/for/results/alphafold_structs \
    --save_directory /path/to/directory/for/results/alphafold_structs/structs \
    --cluster aurum \
    --submit_args ""--job-name=awesome_AF_run"" \
    --skip_existing True
```

# Adding clusters
- All code specific to a computational cluster should be contained within `tps_eval/scripts/<cluster>`.
- Create a file `tps_eval/scripts/<cluster>/config.sh` with cluster-specific settings.
- Create individual job scripts as `tps_eval/scripts/<cluster>/jobs/<job_name>.sh`.

### config.sh
File `tps_eval/scripts/<cluster>/config.sh` needs to include the following:
```sh
SUBMIT_JOB="<command used to submit jobs>" # i.e. "sbatch", "qsub"
```
