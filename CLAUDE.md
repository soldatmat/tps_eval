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
  heavy). **Both structure producers are wired** via `--fold`: `esmfold` (both clusters,
  one whole-FASTA job → `<gen>_esmfold_structs/` + `_pae/`) and `alphafold3` (Aurum-only).
  The AF3 fan-out is an Engine **driver Step** (`Step.driver=True`): it runs
  `scripts/run_alphafold_fanout.sh` *on the login node* (not via sbatch), which calls the
  existing `src/alphafold/run_alphafold_jobs.py` to submit one AF3 job per sequence and
  prints the N job ids; the Engine captures them (`fanout_ids`) so every structure Step
  `afterok`-waits on all N, then an `extract_pae` step populates the PAE dir (PAE-consumers
  wait on it). AF3 holo co-folding via `--af3_cofold {none,mg,mg_ppi,mg_gpp,mg_fpp,mg_ggpp,
  mg_gfpp,mg_ee}`: `mg`/`mg_ppi` place CCD `MG`/`POP`; `mg_<sub>` co-folds one forced prenyl-PP
  substrate (SMILES in `src/alphafold/cofold_substrates.py`) for all designs; `mg_ee` co-folds
  each design's EnzymeExplorer-predicted substrate (the fan-out groups by substrate + a Mg-only
  fallback). With `--fold alphafold3`, mg_ee AUTO-CHAINS EE->cofold: since the login-node fold
  driver can't afterok-wait on the in-pipeline `ee_seq` job, build_steps defers the fold+structure
  branch and submits an `af3_ee_continuation` job (afterok on the whole sequence branch) that
  re-invokes the pipeline with `--enzymeexplorer_csv` once EE is done (scoped via
  `scripts/run_eval_pipeline_continuation.sh`, runs from `args._orig_cwd`, only submits jobs). Pass
  `--enzymeexplorer_csv` only for pre-computed EE / non-in-pipeline folds.
  `scripts/run_alphafold_fanout.sh` builds the input via
  `src/alphafold/build_cofold_input.py` (one CSV per group + manifest) and prints ONE combined
  job-id line. Any non-`none` mode enables the holo tools `ion_site_check` + `substrate_positioning`
  (gated on `run_holo` = cofold!=none OR external structs; `--no_holo_tools` force-skips).
  AF3 ion/ligand placement is a hypothesis — verify at DDXXD/NSE downstream. NOT yet ported (v2):
  EnzymeExplorer-with-structures. Sequence + structure-consuming branches verified end-to-end on
  Aurum; the AF3 fan-out wiring (incl. co-fold modes + holo tools) is dry-run-verified (a live
  fold is expensive); `build_cofold_input` + `substrate_positioning` unit-tested locally.

## Selection & funnel layer (narrowing, NOT metrics)
- `src/selection/` is a SEPARATE layer from the metric tools — it CONSUMES the merged metric table
  and PRODUCES a subset. Its ops are **not** in `pipeline_tools.json` (they aren't metrics) and run
  on a login node (CPU/pandas + mmseqs, no SLURM). Primitives, all CSV-keyed-by-ID:
  `merge.py` (merge per-tool CSVs → wide table, per-cell first-wins + ID-union; mirrors the
  dashboard's merge conventions), `gate.py` (boolean conditions; supports nested `all_of`/`any_of`
  and a `when` clause for class-specific gates), `band_filter.py` (reference bands, per-architecture
  via a `by` categorical), `score.py` (weighted z-sum within group), `diversity_dedup.py` (mmseqs
  best-rep-per-cluster, per-group %id). `select_designs.py` composes them from a JSON spec →
  survivors CSV + FASTA + provenance manifest. `export_bands.py` bridges reference-stats →
  `band_filter` bands_file. Dispatch via `scripts/run_selection.sh <merge|select>`.
  GOTCHA: the composite module is `select_designs.py`, NOT `select.py` — the latter shadows Python's
  stdlib `select` (imported by `subprocess`/`selectors`) and breaks the subprocess calls.
- `scripts/run_funnel.py` chains metric-compute + selection across tiers (stepwise/idempotent; reuses
  `run_eval_pipeline`'s idempotency as the SLURM barrier — see docs/TOOLS.md#run_funnel). Funnel +
  selection recipes are version-controlled under `scripts/funnels/` (`production_300k.json` +
  `select_phase{1,2,3}.json`), verified to reproduce the archived 300k→48 run (Phase-1 exact,
  Phase-2 filter-exact, Phase-3 48/48). Unit tests: `src/selection/test_selection.py`.

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
   Also add the one-liner/branch/default to `scripts/pipeline_tools.json` (powers `--list-tools`).
5. **Document it for end users:** add a `docs/TOOLS.md` section (anchored by tool name —
   purpose/inputs/output columns/method/citation/env+source) and a row in the README "## Tools"
   table linking to that anchor. Keep the README one-liner consistent with the
   `pipeline_tools.json` description. **When the build is split across parallel subagents,
   the agents must NOT edit `README.md`/`docs/TOOLS.md` — the orchestrating session adds all
   the rows/sections in one pass at integration** (concurrent agents clobber these shared files).
6. **Log the change** as an H2 entry in the Obsidian vault `History.md` (see Gotchas —
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
  GOTCHA: the `aggrescan3d` env is an EDITABLE install (`pip install -e vendor/aggrescan3d`);
  a `git submodule update`/reset on that submodule de-registers it (the `aggrescan` launcher
  then dies with `DistributionNotFound`). After any submodule reset, re-run
  `pip install -e vendor/aggrescan3d` in the `aggrescan3d` env.
- **Vendor backup mirrors (insurance against upstream disappearing).** The two
  third-party submodules have independent mirrors (NOT forks) on the `soldatmat` GitHub:
  `github.com/dauparas/ProteinMPNN` → `soldatmat/ProteinMPNN`, and
  `bitbucket.org/lcbio/aggrescan3d` → `soldatmat/aggrescan3d`. (`cif_to_pdb` and
  `pymol_scripts` are already soldatmat-owned.) `.gitmodules` URLs point at upstream —
  the mirrors are pure insurance. They are point-in-time snapshots: **refresh them only
  when you bump a vendored submodule pin** (there's no scheduled sync — nothing else can
  make them stale). One command: `./scripts/refresh_vendor_backups.sh`. If an upstream
  ever dies, swap the URL: `git config -f .gitmodules submodule.vendor/<name>.url
  https://github.com/soldatmat/<name> && git submodule sync && git submodule update --init`.
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
