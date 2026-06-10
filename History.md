# tps_eval — change & metric history

A running log of changes to the pipeline and the evaluation metrics it computes.
Most-recent first. Human-readable companion to `CLAUDE.md` (agent orientation) and
`README.md` (setup/usage). Each metric is a per-design value keyed by `ID`, so they
merge for filtration.

> Convention: append a dated entry whenever a utility/metric is added or changed.
> Keep committed work under "Done"; keep the research-backed backlog under "Planned".

---

## Metric inventory (current)

**Sequence branch**
- `motif_search` — presence of the DDXXD aspartate-rich motif + relaxed acidic
  variants `D[DE]..[DE]`, `[DE][DE]..[DE]`, and the NSE/DTE motif `(N|D)D[LIV].(S|T)...E`.
- `max_sequence_identity` — best % identity to a train/reference set and to self.
- `esm_embedding` + `min_embedding_distance` — ESM-1b embedding, nearest-neighbour
  distance to train/reference and to self.
- `motif_pair_distance` *(new — 2026-06-10)* — residue separation between the
  DDXXD-family and NSE/DTE motifs.

**Structure branch** (consumes a dir of AlphaFold/ESMFold `.pdb`/`.cif`)
- `plddt` — per-structure folding-confidence summary (mean/median/min pLDDT,
  fraction confident, n_residues) read from the B-factor field.
- `structural_identity` *(new — 2026-06-10)* — foldseek best TM-score / LDDT to the
  nearest known-TPS structure (the structural analog of `max_sequence_identity`).
- `motif_structural_distance` *(new — 2026-06-10)* — 3D centroid and min CA-CA
  distance between the two metal-binding motifs' coordinating residues (fold-agnostic).

**Function / expression**
- `enzyme_explorer` (+ `_sequence_only`) — TPS classifier (isTPS + per-product-class scores).
- `soluprot` — solubility prediction.

**Folding (structure producers)**
- AlphaFold3 (Aurum-only, apptainer).
- ESMFold *(in progress — 2026-06-10)*, available on both clusters incl. Karolina.

**Orchestration**
- `scripts/run_eval_pipeline.py` — cluster-agnostic declarative orchestrator
  (supersedes `submit_all.sh` for the sequence + structure-consuming branches).

---

## Done

### 2026-06-10
- **Aggrescan3D aggregation propensity** (`src/structure_metrics/aggregation.py`;
  `vendor/aggrescan3d` submodule) — structure-based aggregation/expressibility signal,
  orthogonal to sequence-based SoluProt. Static mode only (CABS-flex never invoked),
  ~7 s/structure CPU. Output keyed by ID: `a3d_avg_score`, `a3d_total_score`,
  `a3d_max_score`, `a3d_min_score`, `a3d_total_pos_score` (sum of positive = aggregation
  propensity), `n_residues`; `--save_residue_scores` dumps per-residue arrays. Wired into
  the orchestrator structs branch. **Env: Python 2.7** (A3D `master` is Py2-only) +
  bundled freesasa, no FoldX/PyMOL needed — env name `$AGGRESCAN3D_ENV`. Note: AlphaFold
  `OXT` terminal atoms crash A3D's freesasa; the tool strips them via Biopython first.
- **TPS structural-domain composition** (`src/enzyme_explorer/domain_composition.py`) —
  count + types of TPS domains per design via EnzymeExplorer's CPU-only `detect_domains`
  (PyMOL + foldseek, no GPU). Output keyed by ID: `n_domains`, per-type counts for the
  **seven** templates (`alpha, beta, gamma, ids, delta, epsilon, zeta` — `ids` is the
  isoprenyl-diphosphate-synthase-like domain the "six" omits), and `domain_architecture`
  (ordered, e.g. `alpha-beta`). Enumerates the full input ID universe and left-joins, so
  **zero-domain designs still get a row** (`n_domains=0`) — verified on Aurum. Wired into
  the orchestrator structs branch. Uses `$ENZYME_EXPLORER_ENV` (= `enzyme_explorer_prod`).
- **ESMFold structure prediction** (`src/esmfold/`, `scripts/run_esmfold.sh` + per-cluster
  job wrappers; `ESMFOLD_ENV` in `paths.sh`). HuggingFace `transformers`
  `EsmForProteinFolding` (`facebook/esmfold_v1`); writes `<ID>.pdb` mirroring the AF
  structs layout, so `plddt`/`structural_identity`/`motif_structural_distance` consume it
  unchanged — ESMFold confidence needs no new tool. Available on both clusters (unlike
  AF3, Aurum-only). Validated on Aurum (mean pLDDT 77.9, 84% confident on a 380-aa TPS).
  NOTE: HF writes pLDDT on a **0–1 scale**, not 0–100 — `esmfold.py` rescales the
  B-factor ×100 so the 0–100-based `plddt` tool reads it correctly. Env recipe:
  `torch 2.5.1+cu121` + `transformers==4.46.3` (5.x breaks on torch 2.5); no OpenFold.
- **Structure-similarity metrics wired into the orchestrator.**
  - New `structural_identity` tool (`src/structure_metrics/run_structural_identity.py`
    + `scripts/run_structural_identity.sh` + per-cluster job wrappers): foldseek
    best TM-score/LDDT of each generated structure to the nearest known-TPS reference.
    Reference-set-agnostic — point `--known_structs_dir` at MARTS-DB full structures
    or EE reference domains. Smoke-tested on Aurum. *(commit `c78b225`)*
  - Wired `plddt` (`--structs_dir`) and `structural_identity`
    (`--structs_dir --known_structs_dir`) into `run_eval_pipeline.py`; `plots`
    soft-depends on both. *(commit `c78b225`)*
- **Active-site geometry metrics** (`src/structure_metrics/active_site_geometry.py`) —
  catalytic-site-specific, fold-agnostic, apo-robust: `carboxylate_convergence_radius`
  (RMS spread of the DDXXD+NSE/DTE side-chain carboxylate/hydroxyl oxygens about their
  centroid — the putative metal locus), `n_coordinating_oxygens`, `metal_point_void`
  (clearance for the metals), and `catalytic_constellation_rmsd` to reference TPS
  templates (via Biopython Cα+Cβ superposition — PyMOL isn't in the env). Wired into the
  orchestrator structs branch. Discriminates well on Aurum (real 1ps1 synthase 6.75 Å
  vs splayed designs 20–26 Å). fpocket pocket descriptors deferred (installs cleanly).
- **Motif-pair distance metrics (sequence + structure).** Shared motif-localization
  core (`src/sequence_metrics/motif_localization.py`); `motif_pair_distance`
  (sequence) and `motif_structural_distance` (3D, fold-agnostic) tools + orchestrator
  wiring. Smoke-tested on Aurum (205 structures, ~30 Å active-site spans). *(commit `491187a`)*

### Earlier this work cycle (pre-History)
- AlphaFold pLDDT extraction tool (`src/structure_metrics/plddt.py`) + B-factor
  preservation fix in `vendor/cif_to_pdb` (Biopython converter).
- Relaxed DDXXD motif variants + NSE/DTE added to motif search.
- EnzymeExplorer reinstalled from the `revision` branch (`enzyme_explorer_prod`);
  tps_eval EE wrappers re-wired to the new console-script schema; `load_results`
  handles both the old (`isTPS`) and new (`<class>_p_calibrated`) schemas.
- SoluProt setup scripts + `run_soluprot.sh` robustness fixes; plots colour/skip fixes.
- Declarative orchestrator `run_eval_pipeline.py` introduced (single
  `--dependency=afterok:` chain; idempotent).
- Aurum migrated Miniconda → Miniforge; envs relocated.

---

## In progress

- **ESMFold structure prediction** — tool code written (`src/esmfold/`, HF
  `transformers` `EsmForProteinFolding`); writes `<ID>.pdb` with pLDDT in B-factors,
  so the existing `plddt` tool extracts ESMFold confidence unchanged. `esmfold` conda
  env being built + validated on Aurum before commit.
- **TPS structural-domain composition** — extract how many / which TPS domains
  (alpha/beta/gamma/delta/epsilon/zeta) each design has, from EnzymeExplorer's domain
  detection, into a CSV keyed by ID. Edge case handled: designs with **no** detected
  domains still get a row (`n_domains=0`).

---

## Planned metric backlog (research-backed)

Prioritised from a literature/practice survey of de-novo design filters, active-site
metrics, and activity-by-similarity (sources captured in the survey notes). Rank, don't
hard-threshold — absolute cutoffs (e.g. pLDDT) don't transfer across topologies.

**High priority**
- **Self-consistency scRMSD** — ProteinMPNN inverse-fold → refold (ESMFold/AF3) →
  Cα-RMSD back to the design; the standard "designability" metric (accept < 2 Å).
- **ESM pseudo-perplexity (naturalness)** — "One Fell Swoop" single-pass pseudo-PPL on
  the existing ESM stack; plotted against novelty to find the novel-but-protein-like frontier.
- *(done — see Done/2026-06-10)* ~~Active-site carboxylate-cage geometry~~ and
  ~~catalytic-residue constellation RMSD~~.

**Medium priority**
- **Broad "what-else-is-it-like" search** *(in progress)* — DIAMOND/BLAST vs Swiss-Prot
  (sequence) + foldseek vs afdb-swissprot (structure), both annotated so each hit is
  labeled TPS vs non-TPS. Per design: top hit + score, whether it's a TPS, best non-TPS
  hit score, #TPS in top-N. Flags function-drift / confirms specificity. Caveat: TPS
  share folds with prenyltransferases/IDS-like enzymes, so a related-enzyme hit is
  expected — report + annotate, don't auto-penalize.
- **Inter-domain PAE** — relative-orientation confidence between TPS domains (uses the
  EE domain definitions); catches bad two-domain placement pLDDT misses.
- **Active-site pocket descriptors** — fpocket/P2Rank volume/hydrophobicity/enclosure of
  the catalytic cavity vs the natural-TPS distribution (pocket volume tracks product
  class). fpocket installs cleanly in the env; the carboxylate-cage centroid is the
  anchor for selecting the catalytic pocket. *(deferred from the active-site build.)*
- **ProteinMPNN mean NLL** — sequence-given-fold quality (free byproduct of scRMSD).
- **Aromatic / cation-π pocket lining** — count/geometry of Trp/Tyr/Phe stabilizing carbocations.
- **Radius of gyration / compactness** — flag for non-compact predictions (near-free).

**Activity / specificity inference (by similarity to natural TPSs)** — all reuse the
distances tps_eval already computes; emit *coarse* labels with calibrated confidence,
never a high-confidence exact-scaffold prediction.
- **k-NN coarse-class transfer** *(high)* — predict terpene size class (C10/C15/C20/…)
  and cyclic-vs-linear by distance-weighted vote of nearest MARTS-DB neighbours in
  sequence / ESM-embedding / foldseek-structure space, ensembled. Use the annotation-transfer
  cliff as the prior (~40% identity for size class, ≥60% before any scaffold-level claim;
  TM 0.5 fold floor).
- **Specificity-determining-residue (SDR) match + divergence flag** *(high as a negative
  filter)* — match the design's residues at known product-specificity positions to a
  candidate neighbour's signature; raise a flag when a design is globally close to a
  known-product TPS but diverges at SDRs (the TEAS/HPS single-switch regime). The only
  metric that confronts the failure mode head-on.
- **Substrate-class compatibility** *(med-high)* — GPP/FPP/GGPP fit, largely from the
  EnzymeExplorer substrate head + a pocket-volume cross-check.
- **Continuous product-property regression** *(med)* — regress product carbon count
  (strong) and ring count (weak) from embeddings, with prediction intervals.
- **Subfamily / phylogenetic placement** *(med)* — TPS-a/b/c/… clade → coarse class prior;
  also routes which SDR panel the divergence flag uses.

**Low / deferred**
- Substrate docking / AF3 substrate co-folding (wants the Mg/diphosphate context; deep-dive only).
- MolProbity, ddG/ThermoMPNN, instability index — weak discriminators on predicted
  structures or optimization (not triage) tools.

**Caveat carried throughout:** TPS product specificity is *not* robustly predicted by
global similarity — TEAS and HPS are ~80% identical yet make different major products,
switchable by ~4–9 active-site residues. Similarity-transfer metrics are screening aids
for *coarse* properties (is-TPS, size class, substrate, subfamily), not activity guarantees;
fine specificity needs the local SDR/active-site metrics, and even those are weak priors.

---

## Open items
- EE-domains structural-identity run needs domain detection on the generated structures
  first (the domain-composition work above sets this up); the full-structure run against
  MARTS-DB AFDB structures is already supported by `structural_identity`.
- Deploy all changes to Karolina (git pull + submodule update; `paths.sh`
  `ENZYME_EXPLORER_ENV=enzyme_explorer_prod`; install ESMFold env; test orchestrator).
