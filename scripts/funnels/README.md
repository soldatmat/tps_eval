# Funnels — version-controlled selection funnels

A **funnel** narrows a large design pool to a small ordering set through tiers of escalating
compute, applying a selection spec between tiers. `scripts/run_funnel.py` chains the metric
compute (`run_eval_pipeline.py`) and the selection layer (`src/selection/`) tier by tier; the
JSONs here are the version-controlled recipes.

## Files
- `production_300k.json` — **the production 3-phase funnel**: ~325k conditioned-DPLM TPS
  designs → 48 wet-lab candidates (16 per first-cyclization class c0/c1/c10). Tiers:
  1. **phase1** (sequence, no fold) → 2000/class — FPP mono-specificity + isTPS + solubility +
     motifs + novelty gates, ranked by `z(FPP_p_calibrated) − z(sequence_similarity)`.
  2. **phase2** (ESMFold apo) → 100/class — per-architecture pocket bands (by class: c10 =
     single-domain, c0/c1 = two-domain) + shared structure gates + a c10-only geometry gate,
     quality-scored, MMseqs2 diversity-deduped.
  3. **phase3** (AF3 `mg_fpp` holo cofold, Aurum) → 16/class — reference-independent holo gate
     (canonical motif-coordinated Mg + FPP coordination) + novelty gate, quality-ranked
     (`z(mean_plddt)+z(−proteinmpnn_nll)+z(−a3d_avg_score)+z(iptm)`), MMseqs2 diversity-deduped.
  Terminal: order-preparation (yeast codon-opt + Type-3 Golden Gate overhangs).
- `select_phase{1,2,3}.json` — the per-tier selection specs (consumed by `src/selection`).

## Reproduction (acceptance)
Verified against the archived run on NAS `_production_gen_300k_2026-06-21/`:
- **Phase 1** — 309,140 → 6000 (2000/class), **100% ID match** to `survivors_*_rank.csv`.
- **Phase 2** — filter counts **exact** (c0 347 / c1 456 / c10 272 = `*_filtered_ids.csv`);
  diversity-100 overlap 98/81/81 (MMseqs2 clustering picks interchangeable near-siblings).
- **Phase 3** — **48/48 exact** match to `final48.csv` (gate pass c0 91 / c1 96 / c10 88, exact).

Re-run the selection on archived metrics with `run_funnel.py --select-only` (or `run_selection.sh
select --entries … --spec select_phaseN.json`). See `docs/TOOLS.md` → "Selection & funnel".

## Note on cyclization geometry (phase 3)
Cyclization geometry (`foldback_c1_to_distal`, `n_aromatics_lining`) is computed and REPORTED as
a QC sanity check, but is deliberately NOT a filter or ranker in `select_phase3.json`: it is
necessary-not-sufficient with a non-monotonic optimum, and applying it as a hard filter would
actually bind (only 40/48 of the real finalists satisfy `foldback≤7.0 ∧ n_aromatics_lining≥3`),
so it would change the selected set rather than confirm it.
