# Pipeline for in silico evaluation of proteins on computational clusters
- Use `tps_eval/scripts/submit_job.sh` to submit jobs on any cluster.
- Every evaluation tool has its main run script created as `tps_eval/scripts/run_<tool_name>`.
- The main run script can be run from a job script on a computational cluster created as `tps_eval/scripts/<cluster>/jobs/<tool_name>.sh`.

# Installation
You will need to [install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first if you don't have it on your system.

```sh
git clone --recurse-submodules https://github.com/soldatmat/tps_eval
cd tps_eval
./setup.sh
```

If you plan to use SoluProt or EnzymeExplorer calls, redefine the paths to your local installations of the tools and the names of the associated Conda environments in `tps_eval/paths.sh`. You have to install the tools yourself (for SoluProt, see *Optional: SoluProt* below; EnzymeExplorer is installed via its own repo's `scripts/setup_env.sh`).

## Optional: SoluProt
SoluProt (solubility predictor, used by `run_soluprot.sh`) is not a pip/conda package — it's a standalone download plus an old py3.7 conda env and two external binaries (USEARCH, TMHMM). A helper script automates the parts that can be automated:
```sh
./scripts/setup_soluprot.sh [install_dir]        # env + SoluProt code + 64-bit USEARCH
```
It creates the `soluprot` env from `scripts/soluprot_environment.yml` (python 3.7, scikit-learn 0.20.1 — pinned, from the anaconda `defaults` channel), downloads the SoluProt standalone, and fetches a **64-bit** USEARCH v11 (public domain, via `rcedgar/usearch_old_binaries`; the legacy 32-bit `usearch11...i86linux32` will not run on modern x86-64 compute nodes). **TMHMM 2.0 remains manual** (academic license: DTU or `git.loschmidt.cz/misc/tmhmm`) — the script prints the exact unpack location and the `#!/usr/bin/env perl` shebang fix; use its 64-bit `bin/decodeanhmm.Linux_x86_64`. After running, set `SOLUPROT_PATH` / `SOLUPROT_ENV` in `paths.sh`. A Karolina-adapted variant (project-storage paths, cache redirects, `envs_dirs`) lives at `scripts/karolina/setup_soluprot.sh`.

## Optional: switch to licensed PyMOL Incentive
By default `setup.sh` installs `pymol-open-source` (no watermark, no license needed, sufficient for the bundled PyMOL scripts). If you have a Schrödinger PyMOL license and want Incentive features (higher-quality ray-tracing, bundled plugins like APBS), run the add-on script after `setup.sh` to swap `pymol-open-source` for `pymol-bundle`:
```sh
./setup_pymol_bundle.sh                          # install, drop license later
./setup_pymol_bundle.sh --license_path my.lic    # install + auto-place license
```
Without a license, `pymol-bundle` renders a Schrödinger watermark on every image. PyMOL Incentive auto-discovers the license at `$PYMOL_LICENSE_FILE`, `~/.pymol/license.lic`, `<conda_prefix>/share/pymol/license.lic`, or `/etc/pymol/license.lic`; passing `--license_path` just copies your `.lic` into the env's `share/pymol/` so discovery works without env-var setup. A Karolina-specific variant lives at `scripts/karolina/setup_pymol_bundle.sh`.

## macOS specifics
PyMOL currently does not support Native Mac ARM on conda-forge. You will need to create an X86-64 environment and install PyQt5 in addition to the default installation.

```sh
conda create --platform osx-64 -n tps_eval <environment definition from setup.sh>

conda activate tps_eval
# <additional install commands defined in setup.sh>
pip install PyQt5
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
    --no-skip_existing \
    --use_protein_id_as_filename
```

# Adding clusters
- All code specific to a computational cluster should be contained within `tps_eval/scripts/<cluster>`.
- Create a file `tps_eval/scripts/<cluster>/config.sh` with cluster-specific settings.
- Create individual job scripts as `tps_eval/scripts/<cluster>/jobs/<job_name>.sh`.
  - These job scripts should ideally call a cluster-agnostic script `tps_eval/scripts/run_<tool_name>.sh` (or `tps_eval/scripts/<tool_name>/<script_name>.sh` in case there are multiple scripts for running the tool).
- `tps_eval/scripts/submit_job.sh` is a helper script useful for writing cluster-agnostic scripts which take the cluster as an argument and automatically use the correct job submission commands.

### config.sh
File `tps_eval/scripts/<cluster>/config.sh` needs to include the following:
```sh
SUBMIT_JOB="<command used to submit jobs>" # i.e. "sbatch", "qsub"
```
