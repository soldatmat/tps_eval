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
  a needle — amber inside the band, ember outside. With multiple sets the dots carry a
  per-set **outline colour** (pick from a palette in the right panel; a single set has no
  outline). In stratified view the pooled batch rides one "all" row so the per-stratum
  bands stay clean.
- **Manual thresholds**: every numeric column has `min`/`max` inputs (placeholder = the
  p1/p99 default). Set either to override the band used for the in/out colouring; a dashed
  amber guide marks a custom bound.
- **Filter overview** (right panel): per-metric counts of how many designs fall outside
  each column's active threshold (custom or p1–p99), ordered by category like the main
  view; click a row to jump to it. When any custom threshold is set it also shows
  *N / M designs pass all set thresholds*.

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
