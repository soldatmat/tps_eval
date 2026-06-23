# Pipeline for in silico evaluation of proteins on computational clusters
- Use `tps_eval/scripts/submit_job.sh` to submit jobs on any cluster.
- Every evaluation tool has its main run script created as `tps_eval/scripts/run_<tool_name>`.
- The main run script can be run from a job script on a computational cluster created as `tps_eval/scripts/<cluster>/jobs/<tool_name>.sh`.

# Tools
The pipeline is a suite of evaluation tools, each writing a CSV keyed by `ID` so metrics merge for filtration. Run the whole suite with the declarative orchestrator:
```sh
python scripts/run_eval_pipeline.py --cluster <aurum|karolina> --fasta_path gen.fasta \
    [--train_path train.fasta] [--structs_dir structs/] [--known_structs_dir known/]
python scripts/run_eval_pipeline.py --list-tools   # show all tools + select with --only/--include/--exclude
```
The table below summarizes each tool; **full per-tool documentation** (inputs, output columns, method, citations, conda env) is in [`docs/TOOLS.md`](docs/TOOLS.md). Branch is `seq` (keyed off a FASTA → `<input>_<tool>.csv`) or `struct` (keyed off a structures dir → `<structs_dir>_<tool>.csv`).

### Sequence
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [motif_search](docs/TOOLS.md#motif_search) | seq | DDXXD / NSE-DTE motif presence search. | `<fasta>_motifs.csv` |
| [motif_pair_distance](docs/TOOLS.md#motif_pair_distance) | seq | Sequence distance between the two metal-binding motifs. | `<fasta>_motif_pair_distance.csv` |
| [esm_embedding](docs/TOOLS.md#esm_embedding) | seq | ESM-1b embeddings (feeds the min-distance metrics). | `<fasta>_embedding_esm1b.csv` |
| [esm_pseudo_perplexity](docs/TOOLS.md#esm_pseudo_perplexity) | seq | ESM pseudo-perplexity (sequence likelihood / naturalness). | `<fasta>_esm_pseudo_perplexity.csv` |
| [max_sequence_identity](docs/TOOLS.md#max_sequence_identity) | seq | Max pairwise sequence identity (self) / vs the train set. | `<fasta>_max_sequence_identity[_self].csv` |
| [min_embedding_distance](docs/TOOLS.md#min_embedding_distance) | seq | Min ESM-embedding distance (self) / vs train (needs esm). | `<fasta>_embedding_esm1b_min_embedding_distance[_self].csv` |
| [soluprot](docs/TOOLS.md#soluprot) | seq | SoluProt predicted solubility. | `<fasta>_soluprot.csv` |
| [enzyme_explorer_sequence_only](docs/TOOLS.md#enzyme_explorer_sequence_only) | seq | EnzymeExplorer sequence-only TPS classification. | `<fasta>_enzyme_explorer_sequence_only.csv` |
| [swissprot_search](docs/TOOLS.md#swissprot_search) | seq | DIAMOND search vs Swiss-Prot (gen-only; TPS/non-TPS hits). | `<fasta>_swissprot_search.csv` |
| [local_sequence_search](docs/TOOLS.md#local_sequence_search) | seq | Fast LOCAL sequence identity/similarity + top-k neighbours (MMseqs2 default / DIAMOND); feeds the k-NN/SDR sequence space. Complements the global `max_sequence_identity`. | `<fasta>_local_sequence_search.csv` |

### Structure
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [plddt](docs/TOOLS.md#plddt) | struct | AlphaFold/ESMFold pLDDT folding confidence. | `<structs_dir>_plddt.csv` |
| [motif_structural_distance](docs/TOOLS.md#motif_structural_distance) | struct | Structural distance between the two metal-binding motifs. | `<structs_dir>_motif_structural_distance.csv` |
| [active_site_geometry](docs/TOOLS.md#active_site_geometry) | struct | Active-site carboxylate-cage geometry (apo-robust). | `<structs_dir>_active_site_geometry.csv` |
| [radius_of_gyration](docs/TOOLS.md#radius_of_gyration) | struct | Radius of gyration / compactness over Cα atoms. | `<structs_dir>_radius_of_gyration.csv` |
| [pocket_descriptors](docs/TOOLS.md#pocket_descriptors) | struct | Active-site pocket descriptors (fpocket + P2Rank cross-check). | `<structs_dir>_pocket_descriptors.csv` |
| [domain_composition](docs/TOOLS.md#domain_composition) | struct | TPS structural-domain composition (EE CPU detector). | `<structs_dir>_domain_composition.csv` |
| [aggregation](docs/TOOLS.md#aggregation) | struct | Aggrescan3D structure-based aggregation propensity. | `<structs_dir>_aggregation.csv` |
| [foldseek_swissprot_search](docs/TOOLS.md#foldseek_swissprot_search) | struct | Foldseek search vs AlphaFold-Swiss-Prot (TPS/non-TPS hits). | `<structs_dir>_foldseek_swissprot_search.csv` |
| [structural_identity](docs/TOOLS.md#structural_identity) | struct | Foldseek structural identity to nearest known TPS (needs `--known_structs_dir`). | `<structs_dir>_structural_identity.csv` |
| [domain_structural_identity](docs/TOOLS.md#domain_structural_identity) | struct | Domain-level structural identity: EE detects each design's TPS domains, then foldseek-aligns them to the known martsDB reference domains (per-domain-type best TM-score/lddt; `n_detected_domains`). | `<structs_dir>_domain_structural_identity.csv` |
| [proteinmpnn_score](docs/TOOLS.md#proteinmpnn_score) | struct | ProteinMPNN sequence-likelihood (NLL) of the design's own sequence given its fold. | `<structs_dir>_proteinmpnn_score.csv` |
| [self_consistency](docs/TOOLS.md#self_consistency) | struct | HEAVY scRMSD self-consistency (ProteinMPNN → ESMFold refold → RMSD). Opt-in. | `<structs_dir>_self_consistency.csv` |
| [aromatic_lining](docs/TOOLS.md#aromatic_lining) | struct | Aromatic / cation-π pocket lining (Trp/Tyr/Phe count + ring orientation; carbocation-stabilization proxy). | `<structs_dir>_aromatic_lining.csv` |
| [diphosphate_sensor](docs/TOOLS.md#diphosphate_sensor) | struct | Diphosphate-sensor basic residues (Arg/Lys + RY pair) at the metal site. | `<structs_dir>_diphosphate_sensor.csv` |
| [ion_site_check](docs/TOOLS.md#ion_site_check) | struct | Ion-placement check: do AF3 co-folded Mg/Mn ions land in the carboxylate cage? Only carries signal for AF3 holo folds (`--af3_cofold mg*`); apo/ESMFold report `n_ions_modelled=0`. Gated on a holo co-fold / `--no_holo_tools`. | `<structs_dir>_ion_site_check.csv` |
| [substrate_positioning](docs/TOOLS.md#substrate_positioning) | struct | Substrate positioning: is the AF3 co-folded prenyl-PP substrate poised in the cage (diphosphate→cage oxygen/Mg, reactive C1→cage)? Auto-detects the ligand per design (`--af3_cofold mg_<sub>\|mg_ee`); apo/no-substrate → NaN. Gated on a holo co-fold / `--no_holo_tools`. | `<structs_dir>_substrate_positioning.csv` |
| [cyclization_geometry](docs/TOOLS.md#cyclization_geometry) | struct | Cyclization-relevant holo geometry: is the co-folded prenyl-PP substrate folded for cyclization (rgyr, C1→distal fold-back) and lined by an aromatic cation-π track? Reference-independent; necessary-not-sufficient. apo/no-substrate → NaN. Gated on a holo co-fold / `--no_holo_tools`. | `<structs_dir>_cyclization_geometry.csv` |
| [global_confidence](docs/TOOLS.md#global_confidence) | struct | Global fold confidence (pTM/iPTM) from the saved PAE npz (needs `--pae_dir`). | `<structs_dir>_global_confidence.csv` |
| [interdomain_pae](docs/TOOLS.md#interdomain_pae) | struct | Mean/max inter-domain PAE between TPS domains (needs `--pae_dir`; EE domain ranges). | `<structs_dir>_interdomain_pae.csv` |

### Activity / specificity
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [knn_label_transfer](docs/TOOLS.md#knn_label_transfer) | label | Label-agnostic k-NN coarse-label transfer: distance-weighted vote of nearest MARTS-DB neighbours, ensembled across the three similarity spaces, with leave-one-out calibration. Consumes the three `--top_k` CSVs + a `reference_id,label` file. | `<input>_knn_label_transfer.csv` |
| [sdr_divergence](docs/TOOLS.md#sdr_divergence) | struct | Specificity-divergence flag: globally close to a known TPS but divergent at the specificity-determining active-site residues (the TEAS/HPS single-switch regime). | `<structs_dir>_sdr_divergence.csv` |
| [substrate_class](docs/TOOLS.md#substrate_class) | label | Substrate-class combiner: fuses the substrate k-NN vote (3 spaces) with the pocket-volume size band + EnzymeExplorer per-substrate signal into a predicted substrate (GPP/FPP/GGPP/…) + agreement. | `<input>_substrate_class.csv` |

### Folding (structure producers)
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [alphafold3](docs/TOOLS.md#alphafold3) | producer | AlphaFold3 folding (Aurum-only); orchestrator-wired via `--fold alphafold3` (per-sequence fan-out → structs + PAE, then runs the structure branch). | `af_output/` + `structs/<ID>.pdb` + `pae/<ID>_pae.npz` |
| [esmfold](docs/TOOLS.md#esmfold) | producer | ESMFold folding (both clusters); orchestrator-wired via `--fold esmfold` (folds the gen FASTA, then runs the structure branch on the result). | `structs/<ID>.pdb` + `structs_pae/<ID>_pae.npz` |

### Function (structure-dependent)
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [enzyme_explorer](docs/TOOLS.md#enzyme_explorer) | struct | EnzymeExplorer TPS classification with structures; not orchestrator-wired (v2). | `<input>_enzyme_explorer/` |

### Aggregator & visualization
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [plots](docs/TOOLS.md#plots) | aggregator | Merges all enabled metrics into comparison plots. Effectively always on unless excluded. | plot images in `--save_dir` |
| [dashboard](docs/DASHBOARD.md) | aggregator | Builds the interactive, self-contained **natural-bands HTML dashboard** with the design batch overlaid on the committed MARTS-DB reference bands (design metrics with no band are still shown). Default last pipeline step; also standalone via `scripts/run_build_dashboard.sh`. Effectively always on unless excluded. | `dashboard/<gen>_dashboard.html` |
| [plot_domains](docs/TOOLS.md#plot_domains) | visualization | PyMOL images of detected domains overlaid on the structure. Standalone. | PNGs / `.pse` |
| [plot_residue_similarity](docs/TOOLS.md#plot_residue_similarity) | visualization | PyMOL images coloring a design by residue similarity to its matched known structure. Standalone. | PNGs / `.pse` |
| [run_visualization](docs/TOOLS.md#run_visualization) | visualization | Dataset-level class-coloured 2D landscape map of a protein set from any representation (PCA / t-SNE / UMAP / PaCMAP / PCoA; vectors or precomputed distances). Standalone. | PNG figure(s) at `--output` |

### Representation / embedding producers
These emit a per-`id` feature CSV (first column `id`, then feature dims) for a whole protein set; **standalone, NOT orchestrator Steps** — run them directly, then feed into [run_visualization](docs/TOOLS.md#run_visualization).
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [saprot_embedding](docs/TOOLS.md#saprot_embedding) | producer | SaProt-650M structure-aware embeddings (per-residue AA + foldseek-3Di tokens, mean-pooled). Standalone. | feature CSV keyed by `id` |
| [ee_domain_features](docs/TOOLS.md#ee_domain_features) | producer | EnzymeExplorer domain-comparison features (1 − foldseek TM-score to reference functional domains, `domains_subset`). Standalone. | feature CSV keyed by `id` |
| [ee_esm1v_embeddings](docs/TOOLS.md#ee_esm1v_embeddings) | producer | EnzymeExplorer ESM-1v-TPS (subseq-finetuned, 650M) mean embeddings. Standalone. | feature CSV keyed by `id` |
| [active_site_features](docs/TOOLS.md#active_site_features) | producer | 32-d active-site/cation-residue property + geometry profile of the shell within 12 Å of the Mg²⁺ carboxylate cage (class-I; OSC class-12 excluded by label). Standalone. | feature CSV keyed by `id` |

### Reference & orchestration
| Tool | Branch | Description | Output |
|------|--------|-------------|--------|
| [run_eval_pipeline](docs/TOOLS.md#run_eval_pipeline) | orchestrator | Declarative cluster-agnostic orchestrator (config-driven tool selection). | submits SLURM jobs |
| [pipeline_tools.json](docs/TOOLS.md#pipeline_tools) | config | Per-tool default on/off + branch + one-liner driving `--list-tools`. | — |
| [compute_reference_stats](docs/TOOLS.md#compute_reference_stats) | reference | MARTS-DB natural-TPS bands → reference-stats JSON. Standalone. | `src/reference_stats/marts_db_metric_stats.json` |

# Order preparation
A standalone utility (**not** part of the evaluation pipeline) that turns amino-acid
designs into synthesis-ready DNA for Golden Gate / MoClo construction: codon-optimize for
the target organism (default *S. cerevisiae*; internal BsaI/BsmBI removed, homopolymer cap
+ GC window) + add the fixed overhangs. Runs on a login node:
```sh
bash scripts/run_prepare_order.sh designs.fasta        # → designs_order.csv + designs_order.txt
```
Full documentation in [`docs/ORDER_PREPARATION.md`](docs/ORDER_PREPARATION.md).

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
cd "tps_eval/src/alphafold"
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
