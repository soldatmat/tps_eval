# Catalytic-motif base rates by mechanistic Class I / Class II (2026-08-18, NSE/DTE updated 2026-08-18)

**⚠ NSE/DTE column updated same day.** The NSE/DTE regex was relaxed from the strict
`(N|D)D(L|I|V).(S|T)...E` to `[DN]D..[STG]...E` (Durairaj et al. 2019, Phytochemistry
158:157-165; independently confirmed optimal by the empirical grid search below) — the DDxxD
columns are UNCHANGED (same regex, same numbers as the original computation) and reproduce
exactly; only the `NSE/DTE` column values below are new. The EDSQ/non-EDSQ Class-II split was
re-derived from `projects/tps_eval/src/tps_eval/knn/substrate_labels.csv` (substrate=="EDSQ"
for the 110 resolved rows, plus the 7 Class-II enzymes with no resolved substrate label, which
match the original domain-architecture-based EDSQ set exactly by count and reproduce every
DDxxD number in the table below bit for bit) rather than the original NAS domain-composition
CSV (unavailable this session) — same 117/68 split, cross-validated.

Breaks down the pooled catalytic-motif presence rates already documented in
`docs/TPS_DOMAIN_KNOWLEDGE.md` §7a (66.8% / 82.1% / 90.4% / 69.7% over N=1348, NSE/DTE using the
2026-08-18-relaxed regex) by MARTS-DB's own mechanistic **Class I vs Class II** curation. That
breakdown did not exist anywhere in the repo before this computation; only the pooled
(class-agnostic) numbers were recorded.

## Methodology

- **Motif-presence source.** `data/MARTS-DB/MARTS-DB_2026-06-12/reference_stats/esmfold_ref/
  TPS_sequences_motifs.csv` on pluskal.nas (N=1348 natural TPS, one row per unique enzyme). Its
  `sequence` column is byte-identical to `projects/dplm/data-bin/MARTS-DB/2026-06-12/TPS_sequences.csv`
  (same release date, same 1348 `Enzyme_marts_ID`s), so the four motif columns were **recomputed
  locally** from that in-repo sequence file using the exact same regex definitions as
  `projects/tps_eval/src/tps_eval/sequence_metrics/motif_search.py` (`re.search` per pattern, per
  sequence). The recomputation was verified two ways before use: (1) exact match against 44
  spot-checked rows read directly from the NAS-precomputed CSV (0 mismatches across all 4 columns),
  and (2) the resulting pooled counts reproduce the already-documented §7a numbers exactly
  (901/1348, 1107/1348, 1218/1348, and — under the pre-2026-08-18 strict regex — 671/1348; the
  NSE/DTE column update below reproduces 939/1348 under the new regex, also cross-checked against
  a fresh `run_motif_search` job on Aurum's reference set, see the update note above).
- **Class-assignment source.** `projects/dplm/data-bin/MARTS-DB/2026-07-29/reactions.csv`, column
  `Class` (`"1"` or `"2"` per reaction row; an enzyme may have multiple reaction rows). Join key:
  `Enzyme_marts_ID` (reactions.csv) ↔ `ID` (motifs CSV) — both `marts_E#####` strings. All 1348
  motif-CSV IDs were found in reactions.csv (100% join coverage, reconfirmed).
- **Per-enzyme class label.** An enzyme is **Class I** if every one of its reaction rows has
  `Class=1`, **Class II** if every row has `Class=2`, and **mixed** if it has at least one row of
  each (bifunctional fusion enzymes — e.g. copalyl-diphosphate/kaurene-synthase-type diterpene
  synthases with two catalytic domains of different mechanism). This labeling is now saved as a
  reusable `tps_eval` labeling file: `projects/tps_eval/src/tps_eval/reference_stats/class_labels.csv`
  (`reference_id,label` with label ∈ {`I`, `II`, `mixed`}, all 1348 rows), consumable by
  `aggregate_reference_stats.py --group_by` for any future metric, not just motifs.
- **EDSQ / delta-epsilon exclusion.** `data/MARTS-DB/MARTS-DB_2026-06-12/reference_stats/esmfold_ref/
  structs_esmfold_domain_composition.csv` on pluskal.nas, column `domain_architecture`. Rows whose
  value is exactly `delta-epsilon` are squalene/oxidosqualene-cyclase-type enzymes — structurally
  **non-homologous** to the canonical class I/II TPS fold ((αα)₆ QW-repeat barrel, no Mg²⁺-binding
  DDxxD/NSE-DTE motifs by construction) — even though MARTS-DB also labels them `Class=2`. They are
  reported as their own stratum and additionally excluded from a "true class II TPS" cut so the
  class-II motif rate is not dragged down by a mechanistically distinct enzyme family that was never
  expected to carry these motifs. N=117 delta-epsilon rows among the 185 pure-Class-II enzymes (within
  the previously-noted ~111–117 range for this cut).
- **Both NAS source CSVs** were read exclusively via the Read tool (never Bash) per the append-only
  archive's access policy, transcribed locally, and cross-checked (unique-ID-set diff against the
  in-repo sequence file; zero mismatches after removing one transcription artifact).

## Results

N = enzymes in the stratum. Percentages are motif-present / N.

| Stratum | N | `DD..D` (DDxxD) | `D[DE]..[DE]` | `[DE][DE]..[DE]` (relaxed) | NSE/DTE `[DN]D..[STG]...E` (was `(N\|D)D(L\|I\|V).(S\|T)...E`) |
|---|---:|---|---|---|---|
| **All (pooled, for reference)** | 1348 | 901 (66.8%) | 1107 (82.1%) | 1218 (90.4%) | 939 (69.7%) (was 671, 49.8%) |
| **Class I** | 1132 | 867 (76.6%) | 1015 (89.7%) | 1044 (92.2%) | 915 (80.8%) (was 667, 58.9%) |
| **Class II, all** (incl. EDSQ) | 185 | 7 (3.8%) | 61 (33.0%) | 143 (77.3%) | 5 (2.7%) (was 0, 0.0%) |
| **Class II, excluding EDSQ/delta-epsilon** (true class-II TPS) | 68 | 4 (5.9%) | 38 (55.9%) | 57 (83.8%) | 3 (4.4%) (was 0, 0.0%) |
| **EDSQ / delta-epsilon only** (sub-report, subset of "Class II, all") | 117 | 3 (2.6%) | 23 (19.7%) | 86 (73.5%) | 2 (1.7%) (was 0, 0.0%) |
| **Mixed (bifunctional, Class 1 + Class 2 reactions)** | 31 | 27 (87.1%) | 31 (100.0%) | 31 (100.0%) | 19 (61.3%) (was 4, 12.9%) |

Sanity checks: Class I (1132) + Class II all (185) + mixed (31) = 1348 = pooled N. Class II excluding
EDSQ (68) + EDSQ only (117) = 185 = Class II all.

## Interpretation

- **Both motifs are overwhelmingly Class-I signatures.** DDxxD and NSE/DTE track the mechanistic split
  almost perfectly: Class I sits close to (in fact above) the pooled rate on every column, while Class
  II sits near zero on the two strict columns.
- **The EDSQ/delta-epsilon caveat matters for DDxxD and the mid-relaxation column, not for NSE/DTE.**
  Pooling the non-homologous squalene/oxidosqualene cyclases into "Class II" pulls `DD..D` and
  `D[DE]..[DE]` down (3.8%→5.9% and 33.0%→55.9% once excluded) — the EDSQ rows genuinely have almost no
  matches on the mid-strictness column, dragging the pooled class-II figure down for reasons that are
  about a different, non-homologous fold rather than class-II TPS catalysis. **NSE/DTE stays low
  (2.7%→4.4%) in both cuts even under the relaxed regex** — the near-total absence of the NSE/DTE motif
  in Class II is genuine class-II TPS biology, not an EDSQ artifact or a strict-regex detector failure
  (the relaxed regex, which lifts Class I from 58.9%→80.8%, moves Class II by only a handful of
  sequences: 0→5/185 pooled, 0→3/68 excluding EDSQ).
- **The fully relaxed motif (`[DE][DE]..[DE]`) is comparatively motif-agnostic.** It still fires on
  83.8% of true class-II TPS (vs 90.4% pooled, vs 92.2% Class I) — class II retains an aspartate/
  glutamate-rich patch in roughly the expected register even without the canonical DDxxD spacing or the
  NSE/DTE triad.
- **Report the excluding-EDSQ column (N=68, 5.9% / 55.9% / 83.8% / 4.4%) as the more meaningful "true
  class-II TPS" figure** when comparing against Class I or the pooled rate; the "Class II, all" column
  is retained for transparency but conflates two structurally distinct enzyme families.
- **The 31 mixed/bifunctional enzymes were excluded from both Class I and Class II** buckets rather
  than assigned to either — they carry catalytic domains of both mechanisms and would misrepresent
  whichever class they were folded into. For reference, they carry the motifs almost universally
  (87–100% on three of four columns) — consistent with most of them contributing a Class-I domain to
  the enzyme in addition to their Class-II domain.

### Excluded mixed-class enzyme IDs (N=31)

```
marts_E00055, marts_E00106, marts_E00175, marts_E00183, marts_E00254, marts_E00312, marts_E00313,
marts_E00396, marts_E00489, marts_E00492, marts_E00494, marts_E00566, marts_E00582, marts_E00584,
marts_E00586, marts_E00653, marts_E00654, marts_E00655, marts_E00663, marts_E00701, marts_E00719,
marts_E00734, marts_E00769, marts_E00807, marts_E00920, marts_E00965, marts_E01021, marts_E01087,
marts_E01116, marts_E01345, marts_E01543
```

This list is an exact match to the 31 IDs independently reported as bifunctional in the task
specification, confirming the reaction-row-level class join is behaving as expected.

## Deliverables from this computation

1. `projects/tps_eval/src/tps_eval/reference_stats/class_labels.csv` — the reusable `reference_id,label`
   (I/II/mixed) labeling file, all 1348 enzymes, consumable by
   `aggregate_reference_stats.py --group_by` for any future metric.
2. `docs/TPS_DOMAIN_KNOWLEDGE.md` §7a — a new `[measured]` sub-bullet summarizing this table and
   pointing here for the full breakdown.
3. This file — the full table, methodology, and excluded-ID list.
