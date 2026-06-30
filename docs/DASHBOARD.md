# Natural-bands dashboard

An interactive, self-contained HTML view of the MARTS-DB reference bands, with an
optional overlay of a generated-design batch. Standalone reporting utility — not part
of the eval / submit-all pipeline. Runs locally (Python stdlib only, no conda env).

## Run

```sh
# one design set:
scripts/run_build_dashboard.sh --designs '[name=]path[,path2,...]' [--bands <json...>] [--output <out.html>]
# several sets (repeat --designs; each becomes a distinctly-outlined overlay):
scripts/run_build_dashboard.sh --designs 'rfdiffusion=run1/' --designs 'esm3=run2/metrics.csv'
# reference-only / demo:
scripts/run_build_dashboard.sh            # bands only
python3 src/dashboard/build_dashboard.py --demo   # synthetic overlay
```

Each `--designs` value is one set: `[name=]path[,path2,...]`, where each path is a merged
CSV, a directory, or a glob of the pipeline's `*_<tool>.csv` outputs.

## Large batches (10^5+ designs) — "large mode"

Plotting one needle per design dies well before 10^5. **Above 20 000 total designs the builder
auto-switches to large mode** (override with `--large-mode` / `--no-large-mode`, threshold with
`--large-threshold N`). In large mode:

- Per-design values ship as compact **base64 typed arrays** (Float32 / category codes) instead of
  JSON lists, so a single self-contained file stays openable.
- The middle column draws an **aggregate density curve per set + a small sampled scatter** (in/out
  coloured) instead of one needle per design.
- The funnel filters **all** designs via typed arrays (no per-design string set), so the
  *exact* "how many pass" counts stay live while you drag thresholds; counts render compact (`234k`).
- Adjusting one metric's threshold redraws only that card (the funnel still recomputes fully).
- Hide / show / search toggle a CSS class on already-built cards (no DOM rebuild); search does not
  touch the funnel at all. The funnel hot loop uses flattened typed arrays.
- De/selecting a design set (the chips) and recolouring a set redraw only the *visible* cards'
  overlays + the funnel (no card rebuild) — at draw time the bands recompute their axis domain and
  density/needles from the currently-included sets. So the only full rebuilds left are source /
  stratify changes.

Small batches are untouched: plain JSON value lists, per-design needles, classic funnel.

Measured on a 300 000-design synthetic set: ~0.4 s to parse + decode + first render, ~tens of ms
per threshold change. **Payload scales with column count** — large mode keeps every design's value
for every column, so very wide tables (~150 columns × 10^5) produce a large file; ship only the
columns you need to filter on for the biggest batches. Validate quickly with
`python3 src/dashboard/build_dashboard.py --demo --demo-n 300000`.

## Pipeline integration

The dashboard is the **default last step of `run_eval_pipeline.py`** (tool key
`dashboard`, like `plots`): it soft-depends on every enabled metric step and runs once
they finish, overlaying the generated batch on the bands. It writes
`<gen_dir>/dashboard/<gen>_dashboard.html`. Disable with `--exclude dashboard`; it is a
tiny stdlib SLURM job (`scripts/<cluster>/jobs/dashboard.sh`, ~8 GB, 15 min).

**Missing bands are tolerated.** Any design column without a reference band — the
sequence / comparative metrics that have no natural band by design, a freshly-added
metric, or even *all* bands absent — is still shown as a *design-only* card (needles on a
design-derived axis, tagged "no reference band"), so a batch is always inspectable.

- No flags → the reference-only explorer over the committed bands.
- `--demo` → overlays a reproducible synthetic design batch (sanity-checks the overlay).
- `--designs <path>` → overlays a real batch (see below).
- Output defaults to `data/dashboard/marts_db_dashboard.html` (gitignored). Open it in a
  browser; everything (CSS, JS, data) is inlined, so the file is portable and works
  offline / as a published artifact.

## What it shows

- **Sources**: ESMFold / AlphaFold3 / Boltz-2 (holo) — the committed band JSONs in
  `src/reference_stats/marts_db_<source>_metric_stats.json`. Toggle between them.
- **Categories**: metrics are grouped into *Fold & confidence · Active site · Sequence ·
  Function · Novelty* (the comparative similarity metrics), in both the left nav and the
  main view. The grouping lives in `src/dashboard/metric_info.py` (`METRIC_CATEGORY`).
- **Per-metric "?"**: hover the `?` by a metric name for a one-line explanation; numeric
  columns show their **mathematical range** (e.g. `range 0–1`) in the header. Both come
  from `src/dashboard/metric_info.py` (`METRIC_INFO`).
- **Numeric metrics** render as a layered *natural band*: min–max hairline → p1–p99 →
  p5–p95 → p25–p75 (IQR) core → median tick (+ a mean diamond).
- **Categorical metrics** render as stacked-proportion bars of the category frequencies.
- **Stratify by** `substrate` / `first_cyclization` / `domain_architecture` → one small
  band per stratum, sharing the metric's x-axis (top strata by count; the rest summarised).
- **Design overlay** (when ≥1 set loaded): each design's value drops onto the same axis as
  a needle. With a single set the dot is amber inside the band, ember outside; with multiple
  sets each set has its own dot **colour** (full when kept, dimmed-but-coloured when filtered
  out). The **design-set chips** at the top of the Distributions column are the key *and* the
  controls: each shows the set's dot colour + count, a palette to recolour its dots, and is a
  **click-to-include/exclude toggle** — a yellow outline marks an included set (all included
  by default); excluding one drops its designs from every card and from the funnel. The header
  "Designs" box is just a passive total (sets · designs). In stratified view the pooled batch
  rides one "all" row so the per-stratum bands stay clean.
- **Manual thresholds**: every numeric column has `min`/`max` inputs (placeholder = the
  p1/p99 default). Set either to override the band used for the in/out colouring; a dashed
  amber guide marks a custom bound. A small reset button (left of the inputs) clears a custom
  threshold back to the default.
- **Categorical filters**: boolean / categorical columns are filterable too. A
  low-cardinality categorical column (≤ 12 distinct design values) shows a row of **toggle
  pills** — one per category, coloured to match the stacked bar — in place of the min/max
  inputs. All categories start allowed; click a pill to exclude that category (it dims +
  strikes through), and the funnel cuts designs whose value is excluded. The left-hand reset
  re-allows everything. Free-form identifier columns (hit IDs, etc.) exceed the cardinality
  cap and are not offered as filters.
- **Hide / reorder**: hide a metric via the eye toggle in the left column **or** the eye in
  a card's top-right corner (it disappears from both the main view and the pipeline);
  **hide all** hides every metric, **hide no-data** bulk-hides metrics with no design values,
  **show all** restores. Drag a
  card by its grip to reorder the **main** view; drag a pipeline row by its grip to reorder
  a single filter, or drag a pipeline **category header** to move that whole category block
  (this changes the order filters are applied in — see below). Main and pipeline orders are
  independent.
- **Filter pipeline** (right panel) — a filtering funnel, drag-reorderable (click a row to
  jump to its metric). Filters are nested **category › metric subgroup › column**, mirroring
  the left/main columns: each metric's columns sit under a collapsible **subgroup header**
  (click it to collapse the whole metric to a single line showing its aggregate kept/cut;
  drag its grip to move the whole group; click its ✕ to send the whole group to the unused
  tray). Three columns per metric:
  - **PASS** — designs inside that metric's band alone (custom threshold if set, else
    p1–p99); matches the in/out dot colouring.
  - **KEPT** — cumulative survivors applying each metric's band top→bottom (custom
    threshold if set, else p1–p99), with a yellow pipe (its own column between the KEPT and
    CUT numbers) that narrows as designs are filtered out; set/drag a threshold and the
    whole funnel updates live.
  - **CUT** — designs removed at *that step only* (so once the pipe is empty, rows below
    show 0).
  A boxed *designs kept N/M* total sits at the bottom, with an **extract** button that
  downloads the kept designs (set, ID + every metric value) as CSV.
  **Removing a filter from the funnel**: click the **✕** on a row (or drag its grip) to drop
  it into the **unused filters** tray above the funnel; the funnel recomputes without it.
  Re-add it with the chip's **+** button, by dragging the chip back onto the funnel, or
  **use all** to restore every unused filter at once. Removing every filter keeps all designs.

## Overlaying a real design batch

The generator reads the tps_eval pipeline's own per-tool outputs. Point `--designs` at
either a single merged CSV or a **directory of `*_<tool>.csv` files** (the orchestrator's
layout). Every CSV must be keyed by `ID`; columns are matched to the band columns **by
name** (the pipeline already emits the same column names as the bands), so no mapping is
needed — any column that exists in both is overlaid, the rest are ignored.

## Architecture

- `src/dashboard/build_dashboard.py` — loads the band JSON(s) + optional design batch,
  compacts them to the stats the UI needs, and injects the result as inline JSON into
  the template (`/*__DASHBOARD_DATA__*/` token). Non-finite values are sanitised to
  `null` so the embedded JSON is strict.
- `src/dashboard/template.html` — the editable HTML/CSS/JS (vanilla, SVG charts, no
  external dependencies).
- `scripts/run_build_dashboard.sh` — the standard wrapper.
