#!/usr/bin/env python3
"""Stepwise, idempotent driver for a multi-tier selection FUNNEL over the tps_eval pipeline.

A funnel narrows a large design pool to a small ordering set through several tiers of
escalating compute (cheap sequence metrics on everything -> ESMFold on thousands -> AF3 holo
on hundreds), applying a selection spec between tiers. tps_eval already has the per-tier
COMPUTE (run_eval_pipeline.py) and the SELECTION layer (src/selection); this driver is the
thin connective tissue that chains them and carries survivors forward.

It is deliberately STEPWISE (one `--tier N` invocation at a time), mirroring
run_eval_pipeline.py's own idempotent/resumable model and the SLURM reality that a tier's
metrics are asynchronous jobs. The barrier is handled by REUSING run_eval_pipeline's
idempotency: for a tier we run it for real; if it SUBMITS any job, the metrics aren't ready
-> we stop and tell you to monitor + re-run; if it submits NOTHING (every output already
exists), the metrics are ready -> we run the tier's selection and write its survivors (the
next tier's input). No long-lived orchestration process, nothing to leak.

Usage (drive it tier by tier):
    # Tier 0 (sequence metrics), submits jobs then stops:
    python scripts/run_funnel.py --config scripts/funnels/production_300k.json \
        --cluster karolina --workdir RUN/ --seed_fasta all_designs.fasta --tier 0
    # ...monitor the jobs (squeue / the run's logs)...
    # Re-run the SAME command: now outputs exist -> it selects and writes phase1_survivors.*
    # Then tier 1 (esmfold), tier 2 (AF3 holo, on aurum), etc. The last tier's selection is
    # followed by the terminal order-preparation step.

Config schema: see scripts/funnels/production_300k.json. Each tier =
{name, [input], [fold], [af3_cofold], [cluster_override], only:[...tool keys...], select:<spec>.json}.
Selection specs: see scripts/funnels/select_phaseN.json (consumed by src/selection).

Pure stdlib (runs on a login node). It shells out to run_eval_pipeline.py (metrics),
scripts/run_selection.sh (merge + select, in the tps_eval conda env), and
scripts/run_prepare_order.sh (terminal). DO NOT add non-stdlib deps here.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_EVAL = os.path.join(REPO, "scripts", "run_eval_pipeline.py")
RUN_SELECTION = os.path.join(REPO, "scripts", "run_selection.sh")
RUN_PREPARE_ORDER = os.path.join(REPO, "scripts", "run_prepare_order.sh")


def _base(fasta: str) -> str:
    for ext in (".fasta", ".fa"):
        if fasta.endswith(ext):
            return fasta[: -len(ext)]
    return os.path.splitext(fasta)[0]


def metric_globs(fasta: str, fold) -> list:
    """The globs matching a tier's metric CSVs, mirroring run_eval_pipeline's output
    conventions: sequence CSVs sit next to the FASTA (<base>_<tool>.csv); structure CSVs
    sit beside the produced structures dir (<structs>_<tool>.csv)."""
    base = _base(fasta)
    if fold == "esmfold":
        structs = base + "_esmfold_structs"
        return [structs + "_*.csv", base + "_*.csv"]
    if fold == "alphafold3":
        structs = base + "_af3/structs"
        return [structs + "_*.csv", base + "_*.csv"]
    return [base + "_*.csv"]


def tier_input_fasta(cfg: dict, idx: int, args) -> str:
    """Tier 0 reads the seed FASTA; tier k reads the previous tier's survivors FASTA."""
    if idx == 0:
        if not args.seed_fasta:
            raise SystemExit("[FAIL] tier 0 needs --seed_fasta (the design pool to funnel).")
        return os.path.abspath(args.seed_fasta)
    prev = cfg["tiers"][idx - 1]["name"]
    prev_fasta = os.path.join(args.workdir, f"{prev}_survivors.fasta")
    if not os.path.isfile(prev_fasta):
        raise SystemExit(f"[FAIL] tier {idx} needs the previous tier's survivors "
                         f"({prev_fasta}); run --tier {idx - 1} to completion first.")
    return prev_fasta


def run_metrics(tier: dict, fasta: str, cluster: str, args) -> bool:
    """Run run_eval_pipeline for this tier. Returns True if it SUBMITTED any job (metrics
    not yet ready -> caller stops), False if everything was already done (-> select)."""
    cmd = [sys.executable, RUN_EVAL, "--cluster", cluster, "--fasta_path", fasta,
           "--only", ",".join(tier["only"]),
           # The funnel never needs per-tier plots/dashboard; excluding them keeps the
           # "did it submit anything?" signal clean (those aggregators always (re)submit).
           "--exclude", "plots,dashboard"]
    if tier.get("fold"):
        cmd += ["--fold", tier["fold"]]
    if tier.get("af3_cofold"):
        cmd += ["--af3_cofold", tier["af3_cofold"]]
    if args.train_path:
        cmd += ["--train_path", os.path.abspath(args.train_path)]
    if args.known_structs_dir:
        cmd += ["--known_structs_dir", os.path.abspath(args.known_structs_dir)]
    if args.account:
        cmd += ["--account", args.account]
    if args.dry_run:
        cmd += ["--dry-run"]
    print(f"[funnel] metrics: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = (proc.stdout or "") + (proc.stderr or "")
    print(out)
    if proc.returncode != 0:
        raise SystemExit(f"[FAIL] tier metric pipeline failed (rc={proc.returncode}).")
    # A real submission prints '[sub ] ...' or a non-empty '[fan ] ... submitted N ...'.
    submitted = bool(re.search(r"^\[sub \]", out, re.M)) or \
        bool(re.search(r"^\[fan \].*submitted [1-9]", out, re.M))
    if args.dry_run:
        # Dry-run never really submits; treat '[dry ]' lines as "would submit".
        submitted = bool(re.search(r"^\[dry \]", out, re.M))
    return submitted


def run_selection(tier: dict, fasta: str, cfg: dict, args) -> str:
    """Merge the tier's metrics + the previous survivors (carry-forward) and apply the
    tier's selection spec. Returns the survivors CSV path."""
    idx = [t["name"] for t in cfg["tiers"]].index(tier["name"])
    spec = os.path.join(REPO, "scripts", "funnels", tier["select"]) \
        if not os.path.isabs(tier["select"]) else tier["select"]
    entries = metric_globs(fasta, tier.get("fold"))
    if idx > 0:  # carry forward the previous tier's columns (e.g. sequence_identity, plddt)
        prev = cfg["tiers"][idx - 1]["name"]
        prev_csv = os.path.join(args.workdir, f"{prev}_survivors.csv")
        if os.path.isfile(prev_csv):
            entries = [prev_csv] + entries
    out_prefix = os.path.join(args.workdir, tier["name"])
    cmd = ["bash", RUN_SELECTION, "select", "--entries", *entries,
           "--spec", spec, "--output_prefix", out_prefix, "--fasta", fasta,
           "--title", f"{cfg.get('name', 'funnel')} / {tier['name']}"]
    print(f"[funnel] select: {' '.join(cmd)}")
    if args.dry_run:
        print("[funnel] (dry-run) would run selection as above.")
        return out_prefix + "_survivors.csv"
    proc = subprocess.run(cmd, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"[FAIL] tier selection failed (rc={proc.returncode}).")
    return out_prefix + "_survivors.csv"


def run_terminal(cfg: dict, args) -> None:
    term = cfg.get("terminal")
    if not term or term.get("op") != "order_preparation":
        return
    last = cfg["tiers"][-1]["name"]
    survivors_fasta = os.path.join(args.workdir, f"{last}_survivors.fasta")
    if not os.path.isfile(survivors_fasta):
        raise SystemExit(f"[FAIL] terminal order-prep needs {survivors_fasta} "
                         f"(run the last tier to completion first).")
    out_prefix = os.path.join(args.workdir, "order")
    cmd = ["bash", RUN_PREPARE_ORDER, survivors_fasta, "-o", out_prefix,
           "--organism", term.get("organism", "yeast"),
           "--overhang_type", term.get("overhang_type", "Type 3")]
    print(f"[funnel] terminal order-preparation: {' '.join(cmd)}")
    if args.dry_run:
        print("[funnel] (dry-run) would run order-preparation as above.")
        return
    proc = subprocess.run(cmd, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"[FAIL] order-preparation failed (rc={proc.returncode}).")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="Funnel config JSON.")
    p.add_argument("--cluster", required=True, help="Cluster for this tier's metric jobs.")
    p.add_argument("--workdir", required=True, help="Run directory (survivors + manifests).")
    p.add_argument("--tier", type=int, required=True, help="Tier index to advance (0-based).")
    p.add_argument("--seed_fasta", default=None, help="The design pool FASTA (tier 0 only).")
    p.add_argument("--train_path", default=None, help="Reference/train FASTA (forwarded).")
    p.add_argument("--known_structs_dir", default=None, help="Known-TPS structs dir (forwarded).")
    p.add_argument("--account", default=None, help="SLURM account (forwarded to metrics).")
    p.add_argument("--select-only", dest="select_only", action="store_true",
                   help="Skip the metric pipeline; assume outputs exist and just select "
                        "(useful to re-select with a changed spec, or on archived metrics).")
    p.add_argument("--terminal", action="store_true",
                   help="After the last tier is selected, run the terminal order-preparation.")
    p.add_argument("--dry-run", action="store_true", help="Print actions; submit/select nothing.")
    args = p.parse_args()

    args.workdir = os.path.abspath(args.workdir)
    os.makedirs(args.workdir, exist_ok=True)
    with open(args.config) as fh:
        cfg = json.load(fh)
    tiers = cfg["tiers"]
    if not 0 <= args.tier < len(tiers):
        raise SystemExit(f"[FAIL] --tier {args.tier} out of range (0..{len(tiers) - 1}).")
    tier = tiers[args.tier]

    cluster = tier.get("cluster_override", args.cluster)
    if tier.get("cluster_override") and tier["cluster_override"] != args.cluster:
        print(f"[funnel] note: tier '{tier['name']}' config expects cluster "
              f"'{tier['cluster_override']}' (you passed --cluster {args.cluster}); "
              f"using '{cluster}'.")

    survivors_csv = os.path.join(args.workdir, f"{tier['name']}_survivors.csv")
    if os.path.isfile(survivors_csv) and not args.dry_run:
        print(f"[funnel] tier '{tier['name']}' already selected ({survivors_csv}).")
        if args.tier + 1 < len(tiers):
            print(f"[funnel] next: --tier {args.tier + 1}")
        elif args.terminal:
            run_terminal(cfg, args)
        return

    fasta = tier_input_fasta(cfg, args.tier, args)
    print(f"=== funnel '{cfg.get('name', '?')}' tier {args.tier} ('{tier['name']}') "
          f"on {cluster}; input {os.path.relpath(fasta, REPO) if fasta.startswith(REPO) else fasta} ===")

    if not args.select_only:
        if run_metrics(tier, fasta, cluster, args):
            print(f"\n[funnel] tier '{tier['name']}' metrics SUBMITTED. Monitor the jobs, then "
                  f"re-run:\n  python scripts/run_funnel.py --config {args.config} "
                  f"--cluster {args.cluster} --workdir {args.workdir} --tier {args.tier}")
            return
        print(f"[funnel] tier '{tier['name']}' metrics all present -> selecting.")

    run_selection(tier, fasta, cfg, args)
    if args.tier + 1 < len(tiers):
        print(f"[funnel] tier '{tier['name']}' done. Next: --tier {args.tier + 1}")
    else:
        print(f"[funnel] final tier '{tier['name']}' done.")
        if args.terminal:
            run_terminal(cfg, args)


if __name__ == "__main__":
    main()
