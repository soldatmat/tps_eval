# Natural-bands dashboard

An interactive, self-contained HTML view of the MARTS-DB reference bands, with an
optional overlay of a generated-design batch. Standalone reporting utility — not part
of the eval / submit-all pipeline. Runs locally (Python stdlib only, no conda env).

## Run

```sh
scripts/run_build_dashboard.sh [--designs <csv|dir|glob ...>] [--demo] \
    [--bands <json...>] [--output <out.html>]
# or directly:
python3 src/dashboard/build_dashboard.py --demo
```

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

- **Sources**: ESMFold / AlphaFold3 / Boltz-2 (holo) — the three committed band JSONs in
  `src/reference_stats/marts_db_<source>_metric_stats.json`. Toggle between them.
- **Numeric metrics** render as a layered *natural band*: min–max hairline → p1–p99 →
  p5–p95 → p25–p75 (IQR) core → median tick (+ a mean diamond).
- **Categorical metrics** render as stacked-proportion bars of the category frequencies.
- **Stratify by** `substrate` / `first_cyclization` / `domain_architecture` → one small
  band per stratum, sharing the metric's x-axis (top strata by count; the rest are
  summarised, never silently dropped).
- **Design overlay** (when a batch is loaded): each design's raw value drops onto the
  same axis as an amber needle — coral if it falls outside the natural p1–p99 envelope —
  plus an "N/M in p5–p95" chip per column. In stratified view the pooled batch rides one
  "all TPS" row so the per-stratum natural bands stay clean.

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
