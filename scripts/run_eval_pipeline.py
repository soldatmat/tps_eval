#!/usr/bin/env python3
"""Cluster-agnostic, declarative orchestrator for the tps_eval pipeline.

Replaces scripts/<cluster>/submit_all.sh. The pipeline is defined ONCE as a list
of Steps (name, job script, args, output path, dependencies, capability gate);
the engine submits each step's SLURM job in dependency order, skipping any step
whose output already exists (idempotent/resumable) and chaining dependencies as a
single `--dependency=afterok:id1:id2:...` (the old bash emitted multiple
`--dependency=` flags, of which SLURM honored only the last).

Cluster differences live in CLUSTERS (submit command, dependency/log formatting,
capabilities) — to add a cluster, add an entry, not a forked script.

Scope: this v1 covers the SEQUENCE branch + plots (motif, esm, min/max
distance/identity, soluprot, enzyme_explorer_sequence_only). The AlphaFold /
structures / plddt branch is a per-sequence FAN-OUT (run_alphafold_jobs.py
submits one alphafold.sh job per sequence) with structs-dependent EE/plddt steps
hanging off it; that is a planned v2 extension and is NOT ported here yet —
passing --structs_dir warns and runs the sequence branch only.

Usage:
    python scripts/run_eval_pipeline.py --cluster aurum \
        --fasta_path gen.fasta [--train_path train.fasta] \
        [--train_embeddings_path <csv>] [--data_colors dodgerblue goldenrod] [--dry-run]

Pure stdlib — runs on a login node, no conda env needed (it only shells out to
`sbatch` and checks output paths).
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/.. == repo root


# --------------------------------------------------------------------------- #
# Cluster adapters                                                            #
# --------------------------------------------------------------------------- #
def _slurm_dep(job_ids: List[str]) -> List[str]:
    # SLURM ANDs colon-separated job ids in a SINGLE --dependency flag.
    return [f"--dependency=afterok:{':'.join(job_ids)}"] if job_ids else []


def _slurm_log(log_dir: str) -> List[str]:
    return [f"--output={log_dir}/%x.%j.out"]


def _slurm_jobid(submit_stdout: str) -> str:
    m = re.search(r"Submitted batch job (\d+)", submit_stdout)
    if not m:
        raise RuntimeError(f"could not parse a SLURM job id from:\n{submit_stdout}")
    return m.group(1)


CLUSTERS: Dict[str, dict] = {
    "aurum": {"submit": "sbatch", "dep": _slurm_dep, "log": _slurm_log,
              "jobid": _slurm_jobid, "caps": {"alphafold"}},
    "karolina": {"submit": "sbatch", "dep": _slurm_dep, "log": _slurm_log,
                 "jobid": _slurm_jobid, "caps": set()},  # AlphaFold is Aurum-only
}


# --------------------------------------------------------------------------- #
# Output-path conventions (must match the tools' own _get_save_path naming)   #
# --------------------------------------------------------------------------- #
def _base(fasta: str) -> str:
    for ext in (".fasta", ".fa"):
        if fasta.endswith(ext):
            return fasta[: -len(ext)]
    return os.path.splitext(fasta)[0]


def out_motif(f): return _base(f) + "_motifs.csv"
def out_esm(f): return _base(f) + "_embedding_esm1b.csv"
def out_mindist(f): return _base(f) + "_embedding_esm1b_min_embedding_distance.csv"
def out_mindist_self(f): return _base(f) + "_embedding_esm1b_min_embedding_distance_self.csv"
def out_maxid(f): return _base(f) + "_max_sequence_identity.csv"
def out_maxid_self(f): return _base(f) + "_max_sequence_identity_self.csv"
def out_soluprot(f): return _base(f) + "_soluprot.csv"
def out_ee_seq(f): return _base(f) + "_enzyme_explorer_sequence_only.csv"


# --------------------------------------------------------------------------- #
# Step model + engine                                                         #
# --------------------------------------------------------------------------- #
@dataclass
class Step:
    name: str
    job: str                       # job-script basename under scripts/<cluster>/jobs/
    args: List[str]
    output: str                    # file OR dir whose existence means "already done"
    deps: List[str] = field(default_factory=list)        # HARD deps: skip step if unmet
    soft_deps: List[str] = field(default_factory=list)   # wait-on-if-submitted; never gates
    requires_cap: Optional[str] = None


class Engine:
    def __init__(self, cluster: str, log_dir: str, dry_run: bool = False):
        self.cfg = CLUSTERS[cluster]
        self.cluster = cluster
        self.log_dir = log_dir
        self.dry_run = dry_run
        self.jobs_dir = os.path.join(REPO, "scripts", cluster, "jobs")
        self.job_ids: Dict[str, str] = {}   # step -> SLURM id (submitted steps only)
        self.satisfied: set = set()         # steps whose output exists or were submitted

    def _exists(self, path: str) -> bool:
        if os.path.isdir(path):
            return any(os.scandir(path))    # dir outputs count only if non-empty
        return os.path.isfile(path)

    def run(self, steps: List[Step]) -> None:
        for s in steps:
            if s.requires_cap and s.requires_cap not in self.cfg["caps"]:
                print(f"[cap ] {s.name}: skipped ('{s.requires_cap}' unavailable on {self.cluster})")
                continue
            missing = [d for d in s.deps if d not in self.satisfied]
            if missing:
                print(f"[dep ] {s.name}: skipped (unsatisfied deps: {', '.join(missing)})")
                continue
            if self._exists(s.output):
                print(f"[skip] {s.name}: output exists ({os.path.relpath(s.output, REPO)})")
                self.satisfied.add(s.name)
                continue
            job_path = os.path.join(self.jobs_dir, s.job)
            if not os.path.isfile(job_path):
                print(f"[MISS] {s.name}: job script not found ({job_path}) -- skipping")
                continue

            dep_ids = [self.job_ids[d] for d in (s.deps + s.soft_deps) if d in self.job_ids]
            cmd = ([self.cfg["submit"]] + self.cfg["dep"](dep_ids)
                   + self.cfg["log"](self.log_dir) + [job_path] + s.args)
            if self.dry_run:
                wait = f" [after {':'.join(dep_ids)}]" if dep_ids else ""
                print(f"[dry ] {s.name}{wait}: {' '.join(cmd)}")
                self.job_ids[s.name] = f"<{s.name}>"
                self.satisfied.add(s.name)
                continue
            proc = subprocess.run(cmd, capture_output=True, text=True)
            out = (proc.stdout or "") + (proc.stderr or "")
            try:
                jid = self.cfg["jobid"](out)
            except RuntimeError as e:
                raise SystemExit(f"[FAIL] {s.name}: submission failed\n{out}") from e
            dep_note = f" (after {':'.join(dep_ids)})" if dep_ids else ""
            print(f"[sub ] {s.name}: job {jid}{dep_note}")
            self.job_ids[s.name] = jid
            self.satisfied.add(s.name)


# --------------------------------------------------------------------------- #
# Pipeline definition                                                         #
# --------------------------------------------------------------------------- #
def build_steps(args) -> List[Step]:
    gen = args.fasta_path
    train = args.train_path
    steps: List[Step] = []

    if args.structs_dir or args.train_structs_dir:
        print("[warn] --structs_dir/--train_structs_dir given, but the AlphaFold/structures/"
              "plddt branch is not yet ported to run_eval_pipeline.py (it is a per-sequence "
              "fan-out). Running the sequence branch only; use the per-tool job scripts for "
              "AF/structs/plddt for now.")

    datasets = [("gen", gen)] + ([("train", train)] if train else [])

    for tag, fa in datasets:
        is_train = tag == "train"
        steps.append(Step(f"motif_{tag}", "motif_search.sh", ["--fasta_path", fa], out_motif(fa)))
        steps.append(Step(f"esm_{tag}", "esm_embedding.sh", ["--fasta_path", fa], out_esm(fa)))
        steps.append(Step(f"maxid_self_{tag}", "max_sequence_identity.sh",
                          ["--fasta_path", fa] + (["--train"] if is_train else []),
                          out_maxid_self(fa)))
        steps.append(Step(f"mindist_self_{tag}", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(fa)] + (["--train"] if is_train else []),
                          out_mindist_self(fa), deps=[f"esm_{tag}"]))
        steps.append(Step(f"soluprot_{tag}", "soluprot.sh", ["--fasta_path", fa], out_soluprot(fa)))
        steps.append(Step(f"ee_seq_{tag}", "enzyme_explorer_sequence_only.sh",
                          ["--fasta_path", fa], out_ee_seq(fa)))

    if train:
        train_emb = args.train_embeddings_path or out_esm(train)
        steps.append(Step("maxid_gen_vs_train", "max_sequence_identity.sh",
                          ["--fasta_path", gen, "--train_path", train], out_maxid(gen)))
        steps.append(Step("mindist_gen_vs_train", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(gen), "--train_embeddings_path", train_emb],
                          out_mindist(gen), deps=["esm_gen", "esm_train"]))

    # plots: depend SOFTLY on every metric (wait on whatever was submitted; still run
    # if some metrics were skipped/absent — the plots tool skips missing-input targets).
    metric_steps = [s.name for s in steps]
    plot_dir = os.path.join(os.path.dirname(os.path.abspath(gen)), "plots")
    if train:
        plot_args = ["--fasta_paths", train, gen, "--data_names", "train", "generated",
                     "--data_colors", args.data_colors[0], args.data_colors[1],
                     "--save_dir", plot_dir]
    else:
        plot_args = ["--fasta_paths", gen, "--data_names", "generated",
                     "--data_colors", args.data_colors[1], "--save_dir", plot_dir]
    # Sentinel output that never exists -> plots always (re)runs after the metrics.
    steps.append(Step("plots", "plots.sh", plot_args, os.path.join(plot_dir, ".__never__"),
                      soft_deps=metric_steps))
    return steps


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cluster", required=True, choices=sorted(CLUSTERS))
    p.add_argument("--fasta_path", required=True, help="Generated sequences FASTA.")
    p.add_argument("--train_path", default=None, help="Reference/train FASTA (enables comparisons).")
    p.add_argument("--structs_dir", default=None, help="(v2) AlphaFold structures dir.")
    p.add_argument("--train_structs_dir", default=None, help="(v2) train structures dir.")
    p.add_argument("--train_embeddings_path", default=None, help="Precomputed train embeddings CSV.")
    p.add_argument("--data_colors", nargs=2, default=["dodgerblue", "goldenrod"],
                   metavar=("TRAIN", "GEN"), help="Matplotlib colors for train/generated.")
    p.add_argument("--dry-run", action="store_true", help="Print the submission plan; don't submit.")
    args = p.parse_args()

    args.fasta_path = os.path.abspath(args.fasta_path)
    for opt in ("train_path", "structs_dir", "train_structs_dir", "train_embeddings_path"):
        if getattr(args, opt):
            setattr(args, opt, os.path.abspath(getattr(args, opt)))

    ts = datetime.now().strftime("%Y_%m_%d_%H%M%S")
    log_dir = os.path.join(REPO, "logs", f"run_eval_pipeline-{ts}")
    if not args.dry_run:
        os.makedirs(log_dir, exist_ok=True)

    steps = build_steps(args)
    print(f"=== tps_eval pipeline on '{args.cluster}' — {len(steps)} steps; logs: {log_dir} ===")
    Engine(args.cluster, log_dir, dry_run=args.dry_run).run(steps)


if __name__ == "__main__":
    main()
