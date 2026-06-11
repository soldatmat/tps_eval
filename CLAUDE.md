# tps_eval — agent orientation

In-silico evaluation pipeline for proteins on HPC clusters. Human setup/usage is in
`README.md`; this file is the quick map for agentic sessions. Keep it terse and
durable — no cluster *state* (that's per-user), no restatements of the README.

## Architecture — every eval tool follows the same shape
- `scripts/run_<tool>.sh` — cluster-agnostic wrapper: parse args, `cd` to repo root,
  source `paths.sh`, `conda activate "$TPS_EVAL_ENV"`, export the libstdc++ fix, then
  `cd src/<subdir>` and run the python entry.
- `scripts/<cluster>/jobs/<tool>.sh` — thin SLURM wrapper: `cd` to repo root, call
  `run_<tool>.sh "$@"` (the `#SBATCH` header is the only cluster-specific part).
- `src/<subdir>/run_<tool>.py` — argv entry; `src/<subdir>/<tool>.py` — logic.
- **Output is a CSV keyed by `ID`** with the metric column(s), so metrics compose/merge
  for filtration. **Sequence-branch** tools key off the fasta → `<input>_<tool>.csv`;
  **structure-branch** tools key off the structures DIR → `<structs_dir>_<tool>.csv`
  (sibling of the dir). Tools emit RAW numbers only — "natural TPS" bands are computed
  separately by the reference-stats pipeline (see below), not baked into the tool.
- `scripts/submit_job.sh --cluster <c> --job_name <tool> [--job_args ...]` submits the
  right job script. `paths.sh` holds conda env names + external-tool install paths.
- **Full pipeline:** `scripts/run_eval_pipeline.py --cluster <c> --fasta_path <gen>
  [--train_path <train>] [--structs_dir <s>] [--known_structs_dir <k>]
  [--self_consistency]` — the cluster-agnostic *declarative* orchestrator (one Step
  list; idempotent; deps chained as a single `--dependency=afterok:…`). It supersedes
  `scripts/<cluster>/submit_all.sh` and covers the **full sequence branch** + plots and
  the **full structure-consuming branch** (pLDDT, structural-identity, motif-structural-
  distance, active-site geometry, domain composition, aggregation, broad foldseek search,
  ProteinMPNN-NLL, radius-of-gyration; scRMSD is opt-in via `--self_consistency`, it's
  heavy). NOT yet ported (v2): the AlphaFold/ESMFold fan-out that *produces* structures
  and EnzymeExplorer-with-structures — pass `--structs_dir` once structures exist.
  Verified end-to-end on Aurum. Wiring audit: every `run_<tool>.sh` except the producers
  (`esmfold`, EE-with-structures) and the standalone `plot_domains`/`plot_residue_similarity`
  viz is a Step here.

## To add a new metric/tool (the pattern — follow it)
1. `src/<subdir>/<tool>.py` (logic → DataFrame keyed by `ID` → CSV) + `run_<tool>.py` (argv).
   - **Structure-branch tool?** Reuse the canonical loader in `src/structure_metrics/plddt.py`
     (af3-`af_output`-vs-flat-dir auto-detection, ID = filename stem, `<structs_dir>_<tool>.csv`
     naming) — `motif_structural_distance.py`/`active_site_geometry.py`/`aggregation.py` all mirror it.
   - **Need the DDXXD / NSE/DTE motif positions?** Use the shared
     `src/sequence_metrics/motif_localization.py` (the single source of truth — regexes +
     coordinating-residue offsets). Don't re-encode the motifs.
2. `scripts/run_<tool>.sh` — copy an existing one (e.g. `run_max_sequence_identity.sh`):
   keep the `conda activate "$<ENV>"` + `export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:..."` block.
   Pick the right env (see Gotchas — not every tool uses `TPS_EVAL_ENV`).
3. `scripts/<cluster>/jobs/<tool>.sh` per cluster — Aurum uses `--constraint=gen-a` for CPU,
   `--constraint=gen-b --gres=gpu:geforce_rtx_3090:1` for GPU (NО `-p`); Karolina uses
   `--partition=qcpu`/`qgpu`. `--time` + `--mem` mandatory.
4. **Wire it into `scripts/run_eval_pipeline.py`**: add an `out_<tool>()` helper + a `Step(...)`
   in `build_steps` (sequence-branch in the per-dataset loop; structure-branch in the
   `if structs:` block). Dry-run to confirm it appears, then add its column(s) to the plots
   (`src/plot/constants.py`) and the reference-stats pipeline if it's an intrinsic property.
5. **Log the change** as an H2 entry in the Obsidian vault `History.md` (see Gotchas —
   NOT a repo file), and commit the code per-tool.

## Gotchas (not visible from the code)
- **`paths.sh` is per-install and is NEVER committed with cluster-specific paths.** Each
  cluster edits it locally (env names, `SOLUPROT_PATH`, `ENZYME_EXPLORER_PATH`).
- The `LD_LIBRARY_PATH="$CONDA_PREFIX/lib:..."` line in the runners is the Karolina
  compute-node libstdc++/`GLIBCXX_3.4.29` fix — keep it in every env-activating runner.
- **AlphaFold pLDDT lives in the B-factor field.** AF3's authoritative source is
  `af_output/<job>/<job>_model.cif` (`B_iso`). `run_plddt` auto-detects an `af_output`
  dir (reads the top-ranked model, keyed by job name) or a flat dir of structures whose
  B-factor already holds pLDDT. The extracted `structs/*.pdb` carry pLDDT only because
  `vendor/cif_to_pdb` was patched (Biopython, preserves B-factor) — structs extracted by
  the old Open Babel converter are zeroed, so prefer `af_output` when unsure.
- **Multiple conda envs — not every tool uses `TPS_EVAL_ENV`.** `paths.sh` names them:
  `TPS_EVAL_ENV` (most tools + foldseek + DIAMOND + ESM), `ESMFOLD_ENV` (ESMFold *and*
  ProteinMPNN — `PROTEINMPNN_ENV` defaults to it), `AGGRESCAN3D_ENV` (**Python 2.7** — A3D
  upstream is Py2-only), `ENZYME_EXPLORER_ENV` (= `enzyme_explorer_prod` on Aurum, the
  `revision` branch), `SOLUPROT_ENV`. Each `run_<tool>.sh` activates the right one.
- **SoluProt / EnzymeExplorer are external installs**, not pip/conda packages. SoluProt:
  `scripts/setup_soluprot.sh` (+ external USEARCH/TMHMM, see README "Optional: SoluProt").
  EnzymeExplorer: its own repo's `scripts/setup_env.sh`. tps_eval only calls them via the
  paths/env names in `paths.sh`.
- `vendor/` holds git submodules — `cif_to_pdb`, `pymol_scripts`, `aggrescan3d` (Py2.7),
  `ProteinMPNN` (ships its own weights): `git submodule update --init --recursive`.
- **`/data/` is gitignored.** Committable reference artifacts therefore live under `src/`,
  NOT `data/` — e.g. `src/homology_search/tps_uniprot_accessions.txt` (the TPS-accession
  classification set), and the reference-stats JSON. Large DBs (Swiss-Prot/afdb-swissprot,
  AFDB structures) live OUTSIDE the repo on each cluster, pointed to by `paths.sh`.
- **Folding:** AlphaFold3 is **Aurum-only**; ESMFold (`run_esmfold.sh`) runs on both
  clusters. ESMFold writes pLDDT on a **0–1 scale** — `esmfold.py` rescales the B-factor
  ×100 so the 0–100-based `plddt` tool reads it. Both produce a `structs/` dir of
  `<ID>.pdb` consumed unchanged by the structure branch.
- **Aurum GPU routing:** use `--constraint=gen-b --gres=gpu:geforce_rtx_3090:1`. `gen-a`+`gpu:1`
  routes to the single-node `a36_96_gpu` partition (`a233`), which is frequently down →
  jobs stuck PENDING (this bit `esm_pseudo_perplexity`). See the `aurum-connect` skill.
- **Project history/decisions are NOT in this repo.** They live in the user's Obsidian
  vault at `/Users/soldatmat/Documents/notes/terpene_generation/History.md` (and `Runs.md`,
  `Ideas.md`), maintained via the `obsidian-history-notes` / `obsidian-notes` workflow.
  Do NOT create a `History.md`/`Runs.md` changelog inside this code repo — record
  project-level changes/decisions as an H2 entry in the vault instead.
