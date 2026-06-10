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
  for filtration. Save name follows `<input>_<tool>.csv`.
- `scripts/submit_job.sh --cluster <c> --job_name <tool> [--job_args ...]` submits the
  right job script. `paths.sh` holds conda env names + external-tool install paths.
- **Full pipeline:** `scripts/run_eval_pipeline.py --cluster <c> --fasta_path <gen>
  [--train_path <train>]` — the cluster-agnostic *declarative* orchestrator (one Step
  list; idempotent; deps chained as a single `--dependency=afterok:…`). Covers the
  sequence branch + plots + the structure-consuming metrics (pLDDT via `--structs_dir`,
  foldseek structural-identity via `--structs_dir --known_structs_dir`); it supersedes
  `scripts/<cluster>/submit_all.sh`. NOT yet ported (v2): the AlphaFold fan-out that
  *produces* structures and EnzymeExplorer-with-structures — `submit_all.sh` still
  carries the fold step; pass `--structs_dir` to the orchestrator once structures exist.

## To add a new metric/tool (the pattern — follow it)
1. `src/<subdir>/<tool>.py` (logic → DataFrame keyed by `ID` → CSV) + `run_<tool>.py` (argv).
2. `scripts/run_<tool>.sh` — copy an existing one (e.g. `run_max_sequence_identity.sh`):
   keep the `conda activate` + `export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:..."` block.
3. `scripts/<cluster>/jobs/<tool>.sh` per cluster — Aurum uses `--constraint=gen-a`
   (NО `-p`); Karolina uses `--partition=qcpu`/`qgpu`. `--time` + `--mem` are mandatory.

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
- **SoluProt / EnzymeExplorer are external installs**, not pip/conda packages. SoluProt:
  `scripts/setup_soluprot.sh` (+ external USEARCH/TMHMM, see README "Optional: SoluProt").
  EnzymeExplorer: its own repo's `scripts/setup_env.sh`. tps_eval only calls them via the
  paths/env names in `paths.sh`.
- `vendor/` holds git submodules (`cif_to_pdb`, `pymol_scripts`):
  `git submodule update --init --recursive`.
- AlphaFold jobs are currently configured for the IOCB **Aurum** cluster only.
