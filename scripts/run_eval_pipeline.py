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

Scope: covers the SEQUENCE branch + plots (motif, esm, min/max distance/identity,
soluprot, enzyme_explorer_sequence_only) AND the structure-level metrics that
consume already-folded structures: pLDDT (folding confidence) and foldseek
structural identity to the nearest known TPS (--structs_dir / --known_structs_dir).
What is NOT yet ported (v2): the AlphaFold per-sequence FAN-OUT that *produces*
those structures, and EnzymeExplorer-with-structures — submit_all.sh still carries
the fold step; pass --structs_dir here once structures exist.

Usage:
    python scripts/run_eval_pipeline.py --cluster aurum \
        --fasta_path gen.fasta [--train_path train.fasta] \
        [--structs_dir structs/] [--known_structs_dir known_tps_structs/] \
        [--train_embeddings_path <csv>] [--data_colors dodgerblue goldenrod] [--dry-run]

Config-driven tool selection
-----------------------------
Which metric tools run is decided by a tools config — `scripts/pipeline_tools.json`
(loaded via stdlib `json`; if absent/unreadable, the built-in DEFAULT_TOOLS table is
used). It maps each tool KEY -> {default: bool, branch: "sequence"|"structure",
description: str}. A KEY may back several Steps (e.g. `motif` -> motif_gen +
motif_train); enabling/disabling the key affects all of them. Each Step is tagged
with its `tool` key so the engine can filter. `self_consistency` defaults OFF (heavy
scRMSD); `plots` is the aggregator and is effectively always on unless excluded.

CLI overrides (applied over the config defaults, in this precedence order):
    --only A,B,...      run ONLY these tool keys (+ plots unless also excluded).
    --include A,B,...   force-enable on top of defaults (e.g. --include self_consistency).
    --exclude A,B,...   force-disable.
    --list-tools        print every tool key (default/branch/description) and exit.
Unknown tool keys are rejected with the valid list. `--self_consistency` is kept as a
back-compat alias for `--include self_consistency`.

To register a NEW tool (e.g. pocket_descriptors): add one entry to
scripts/pipeline_tools.json (and mirror it in DEFAULT_TOOLS below for the
config-missing fallback), then add its `out_<tool>()` helper + a Step tagged
`tool="<key>"` in build_steps. The engine/CLI pick it up with no further changes.

Pure stdlib — runs on a login node, no conda env needed (it only shells out to
`sbatch` and checks output paths). DO NOT add a non-stdlib dependency (no pyyaml);
the config is JSON on purpose.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/.. == repo root
TOOLS_CONFIG = os.path.join(REPO, "scripts", "pipeline_tools.json")

# Single pipeline-wide top-k for the three feeder tools (local_sequence_search,
# min_embedding_distance, structural_identity). Their *_topk.csv outputs feed the
# k-NN label-transfer and SDR-divergence consumers. Built-in fallback; overridden by
# the JSON config "_settings.top_k" and then by the --top_k CLI flag.
DEFAULT_TOP_K = 15

# Default k-NN reference artifacts (committed under src/, see CLAUDE.md "/data/ is
# gitignored" -> committable reference artifacts live under src/).
DEFAULT_KNN_LABEL_FILE = os.path.join(REPO, "src", "knn", "first_cyclization_labels.csv")
DEFAULT_KNN_CALIBRATION = os.path.join(
    REPO, "src", "reference_stats", "knn_calibration_first_cyclization.json")


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
# top-k feeder outputs (consumed by knn / sdr_divergence). The feeder tools key
# the top-k CSV off the SAME input as their single-best CSV:
#   max_sequence_identity  -> <fasta>_max_sequence_identity_topk.csv
#   min_embedding_distance -> <embeddings_csv>_min_embedding_distance_topk.csv
#   structural_identity    -> <structs_dir>_structural_identity_topk.csv
def out_maxid_topk(f): return _base(f) + "_max_sequence_identity_topk.csv"
def out_mindist_topk(emb_csv): return os.path.splitext(emb_csv)[0] + "_min_embedding_distance_topk.csv"
# local_sequence_search (MMseqs2/DIAMOND): the FAST local counterpart of
# max_sequence_identity. Keyed off the fasta. Its _topk.csv is the SEQUENCE-space
# feeder for knn + sdr_divergence (replacing the slow Biopython max_sequence_identity
# for that role; max_sequence_identity stays as the plain GLOBAL novelty metric).
def out_local_search(f): return _base(f) + "_local_sequence_search.csv"
def out_local_search_topk(f): return _base(f) + "_local_sequence_search_topk.csv"
# Self mode (within-set novelty) writes to a distinct _self path (passed explicitly via
# --save_path) so it does not collide with the gen-vs-train output, which shares the
# default <input>_local_sequence_search.csv name. Mirrors maxid_self vs maxid.
def out_local_search_self(f): return _base(f) + "_local_sequence_search_self.csv"
def out_soluprot(f): return _base(f) + "_soluprot.csv"
def out_ee_seq(f): return _base(f) + "_enzyme_explorer_sequence_only.csv"
def out_motif_pair(f): return _base(f) + "_motif_pair_distance.csv"
def out_swissprot_search(f): return _base(f) + "_swissprot_search.csv"
def out_esm_ppl(f): return _base(f) + "_esm_pseudo_perplexity.csv"
# Structure-branch outputs are keyed by the structures DIRECTORY, not the fasta:
# the tools save "<structs_dir>_<tool>.csv" as a sibling of the directory.
def out_plddt(d): return d.rstrip(os.sep) + "_plddt.csv"
def out_structural_identity(d): return d.rstrip(os.sep) + "_structural_identity.csv"
def out_structural_identity_topk(d): return d.rstrip(os.sep) + "_structural_identity_topk.csv"
# Consumer (dependency-heavy) outputs keyed by the gen fasta / structs dir:
def out_knn(f): return _base(f) + "_knn_label_transfer.csv"
def out_sdr_divergence(d): return d.rstrip(os.sep) + "_sdr_divergence.csv"
def out_motif_struct(d): return d.rstrip(os.sep) + "_motif_structural_distance.csv"
def out_active_site_geom(d): return d.rstrip(os.sep) + "_active_site_geometry.csv"
def out_radius_of_gyration(d): return d.rstrip(os.sep) + "_radius_of_gyration.csv"
def out_pocket_descriptors(d): return d.rstrip(os.sep) + "_pocket_descriptors.csv"
def out_domain_composition(d): return d.rstrip(os.sep) + "_domain_composition.csv"
def out_aggregation(d): return d.rstrip(os.sep) + "_aggregation.csv"
def out_foldseek_swissprot(d): return d.rstrip(os.sep) + "_foldseek_swissprot_search.csv"
def out_proteinmpnn(d): return d.rstrip(os.sep) + "_proteinmpnn_score.csv"
def out_self_consistency(d): return d.rstrip(os.sep) + "_self_consistency.csv"
def out_interdomain_pae(d): return d.rstrip(os.sep) + "_interdomain_pae.csv"
def out_global_confidence(d): return d.rstrip(os.sep) + "_global_confidence.csv"
def out_aromatic_lining(d): return d.rstrip(os.sep) + "_aromatic_lining.csv"
def out_diphosphate_sensor(d): return d.rstrip(os.sep) + "_diphosphate_sensor.csv"


# --------------------------------------------------------------------------- #
# Tool catalog (config-driven on/off)                                         #
# --------------------------------------------------------------------------- #
# Built-in fallback used when scripts/pipeline_tools.json is missing/unreadable.
# Keep in sync with that file. To register a NEW tool, add ONE entry here (and to
# the JSON), then tag its Step(s) with tool="<key>" in build_steps -- e.g.:
#   "pocket_descriptors": {"default": True, "branch": "structure",
#                          "description": "Pocket descriptors ..."},
DEFAULT_TOOLS: Dict[str, dict] = {
    "motif":                {"default": True,  "branch": "sequence",  "description": "DDXXD / NSE-DTE motif presence search."},
    "motif_pair":           {"default": True,  "branch": "sequence",  "description": "Sequence distance between the two metal-binding motifs."},
    "esm":                  {"default": True,  "branch": "sequence",  "description": "ESM-1b embeddings (feeds mindist_*)."},
    "esm_ppl":              {"default": True,  "branch": "sequence",  "description": "ESM pseudo-perplexity (sequence likelihood)."},
    "maxid_self":           {"default": True,  "branch": "sequence",  "description": "Max pairwise sequence identity within the dataset."},
    "local_sequence_search": {"default": True, "branch": "sequence",  "description": "Fast LOCAL (MMseqs2) sequence identity/similarity/coverage search; self (within-set novelty) + gen-vs-train. Its _topk.csv is the k-NN/SDR sequence-space feeder."},
    "mindist_self":         {"default": True,  "branch": "sequence",  "description": "Min ESM-embedding distance within the dataset (needs esm)."},
    "soluprot":             {"default": True,  "branch": "sequence",  "description": "SoluProt predicted solubility."},
    "ee_seq":               {"default": True,  "branch": "sequence",  "description": "EnzymeExplorer sequence-only TPS classification."},
    "swissprot_search":     {"default": True,  "branch": "sequence",  "description": "DIAMOND search vs Swiss-Prot (gen-only; TPS/non-TPS hits)."},
    "maxid_gen_vs_train":   {"default": True,  "branch": "sequence",  "description": "Max sequence identity of each gen seq vs the train set."},
    "mindist_gen_vs_train": {"default": True,  "branch": "sequence",  "description": "Min ESM-embedding distance of gen vs train (needs esm)."},
    "plddt":                {"default": True,  "branch": "structure", "description": "AlphaFold/ESMFold pLDDT folding confidence."},
    "motif_struct":         {"default": True,  "branch": "structure", "description": "Structural distance between the two metal-binding motifs."},
    "active_site_geom":     {"default": True,  "branch": "structure", "description": "Active-site carboxylate-cage geometry (apo-robust)."},
    "aromatic_lining":      {"default": True,  "branch": "structure", "description": "Aromatic / cation-pi pocket lining (carbocation-stabilization proxy)."},
    "diphosphate_sensor":   {"default": True,  "branch": "structure", "description": "Diphosphate-sensor basic residues (Arg/Lys + RY pair) at the metal site."},
    "radius_of_gyration":   {"default": True,  "branch": "structure", "description": "Radius of gyration / compactness over Ca atoms."},
    "pocket_descriptors":   {"default": True,  "branch": "structure", "description": "Active-site pocket descriptors (fpocket volume/hydrophobicity/enclosure + P2Rank ligandability cross-check)."},
    "domain_composition":   {"default": True,  "branch": "structure", "description": "TPS structural-domain composition (EE CPU detector)."},
    "aggregation":          {"default": True,  "branch": "structure", "description": "Aggrescan3D structure-based aggregation propensity."},
    "foldseek_swissprot":   {"default": True,  "branch": "structure", "description": "Foldseek search vs AlphaFold-Swiss-Prot (TPS/non-TPS hits)."},
    "structural_identity":  {"default": True,  "branch": "structure", "description": "Foldseek structural identity to nearest known TPS (needs --known_structs_dir)."},
    "proteinmpnn":          {"default": True,  "branch": "structure", "description": "ProteinMPNN sequence-likelihood (NLL) of the design's own sequence given its fold."},
    "global_confidence":    {"default": True,  "branch": "structure", "description": "Global fold confidence (pTM/iPTM) from the saved PAE npz (needs --pae_dir)."},
    "interdomain_pae":      {"default": True,  "branch": "structure", "description": "Mean/max inter-domain PAE between TPS domains (needs --pae_dir; EE domain ranges)."},
    "self_consistency":     {"default": False, "branch": "structure", "description": "HEAVY scRMSD self-consistency (ProteinMPNN -> ESMFold refold -> RMSD). Opt-in."},
    "knn":                  {"default": True,  "branch": "structure", "description": "k-NN coarse-label transfer: ensembled vote over the sequence/embedding/structural top-k neighbours -> predicted coarse label + confidence (needs --train_path + --structs_dir + --known_structs_dir)."},
    "sdr_divergence":       {"default": True,  "branch": "structure", "description": "SDR specificity-divergence: flags designs globally close to a known TPS but divergent at the specificity-determining active-site residues (needs --structs_dir + --known_structs_dir)."},
    "plots":                {"default": True,  "branch": "sequence",  "description": "Aggregator: merges all enabled metrics into plots. Effectively always on unless excluded."},
}


def load_tools_catalog() -> Dict[str, dict]:
    """Load the tools catalog from JSON, falling back to DEFAULT_TOOLS.

    JSON keys that exist in DEFAULT_TOOLS override the built-in entry; any
    DEFAULT_TOOLS key missing from the JSON is filled from the fallback (so a stale
    config that simply lacks a freshly-added tool still runs that tool). Keys
    starting with '_' (e.g. '_comment') are ignored.
    """
    catalog = {k: dict(v) for k, v in DEFAULT_TOOLS.items()}
    try:
        with open(TOOLS_CONFIG) as fh:
            data = json.load(fh)
    except (OSError, ValueError) as e:
        print(f"[note] tools config not used ({TOOLS_CONFIG}: {e}); using built-in defaults.")
        return catalog
    for key, spec in data.items():
        if key.startswith("_"):
            continue
        if not isinstance(spec, dict):
            continue
        entry = catalog.get(key, {"default": True, "branch": "sequence", "description": ""})
        entry.update({k: spec[k] for k in ("default", "branch", "description") if k in spec})
        catalog[key] = entry
    return catalog


def knn_calibrated_spaces(calibration_path: str) -> set:
    """Return the similarity spaces the k-NN calibration JSON actually covers.

    `transfer_labels` indexes calibration["spaces"][space] for EVERY space whose
    top-k CSV is passed, so passing an UNCALIBRATED space (e.g. 'sequence' while only
    embedding+structural are calibrated) raises KeyError. The orchestrator therefore
    only forwards the topk flags for the calibrated spaces. Returns an empty set if the
    file is unreadable (the caller then [note]-skips knn).
    """
    try:
        with open(calibration_path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return set()
    spaces = data.get("spaces")
    return set(spaces) if isinstance(spaces, dict) else set()


def load_top_k() -> int:
    """Resolve the pipeline-wide top-k from JSON "_settings.top_k", else DEFAULT_TOP_K.

    The tools-catalog loader ignores '_'-prefixed keys, so the top-k setting is read
    here explicitly from "_settings": {"top_k": N}. Missing/unreadable/malformed ->
    the built-in DEFAULT_TOP_K fallback.
    """
    try:
        with open(TOOLS_CONFIG) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return DEFAULT_TOP_K
    settings = data.get("_settings")
    if isinstance(settings, dict):
        val = settings.get("top_k")
        if isinstance(val, int) and val >= 1:
            return val
    return DEFAULT_TOP_K


def _parse_keys(csv: Optional[str]) -> List[str]:
    return [t.strip() for t in csv.split(",") if t.strip()] if csv else []


def resolve_enabled_tools(catalog: Dict[str, dict], args) -> set:
    """Compute the set of ENABLED tool keys from config defaults + CLI overrides.

    Precedence: config defaults -> --only (restrict) -> --include (add) -> --exclude
    (remove). `plots` is auto-kept on unless it is explicitly excluded (or --only is
    given without it AND it is not later excluded -> we re-add it).
    Unknown tool keys raise SystemExit listing the valid keys.
    """
    valid = set(catalog)
    only = _parse_keys(args.only)
    include = _parse_keys(args.include)
    exclude = _parse_keys(args.exclude)
    if args.self_consistency:                       # back-compat alias
        include.append("self_consistency")

    unknown = [k for k in (only + include + exclude) if k not in valid]
    if unknown:
        raise SystemExit(
            f"[FAIL] unknown tool key(s): {', '.join(sorted(set(unknown)))}\n"
            f"       valid keys: {', '.join(sorted(valid))}")

    if only:
        enabled = set(only)
        enabled.add("plots")                        # aggregator stays on unless excluded
    else:
        enabled = {k for k, v in catalog.items() if v.get("default", False)}
    enabled |= set(include)
    enabled -= set(exclude)
    return enabled


# --------------------------------------------------------------------------- #
# Step model + engine                                                         #
# --------------------------------------------------------------------------- #
@dataclass
class Step:
    name: str
    job: str                       # job-script basename under scripts/<cluster>/jobs/
    args: List[str]
    output: str                    # file OR dir whose existence means "already done"
    tool: str = ""                 # tool KEY this step belongs to (config-driven on/off)
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
def build_steps(args, enabled: set) -> List[Step]:
    """Build the full candidate Step list, then keep only steps whose tool key is
    enabled. `plots` (the aggregator) soft-depends only on the metric steps that
    survived filtering, so it never waits on a disabled/absent tool. Hard-dep
    semantics are unchanged: if a step's enabled but its hard dep (e.g. esm for
    mindist_*) was filtered out, the Engine emits the existing '[dep ]' skip
    message rather than auto-pulling the dep -- see resolve_enabled_tools docs."""
    gen = args.fasta_path
    train = args.train_path
    steps: List[Step] = []

    structs = args.structs_dir
    known_structs = args.known_structs_dir
    pae_dir = args.pae_dir
    if args.train_structs_dir:
        print("[warn] --train_structs_dir is not used by the orchestrator yet "
              "(EnzymeExplorer-with-structures and the AlphaFold fan-out are not ported).")

    datasets = [("gen", gen)] + ([("train", train)] if train else [])

    # The k-NN + SDR consumers read the gen-vs-train / gen-vs-known *_topk.csv CSVs
    # emitted by the three feeders: local_sequence_search (sequence, MMseqs2 — fast),
    # min_embedding_distance (embedding), structural_identity (structural). Only ask the
    # feeders for top-k when a consumer is actually enabled (no wasted top-k otherwise).
    # topk_args is appended to the three gen-vs-{train,known} feeder Steps below.
    # max_sequence_identity is NO LONGER a feeder — it stays as the plain GLOBAL novelty
    # metric (no --top_k); the fast local search supplies the sequence-space top-k.
    topk_consumers = {"knn", "sdr_divergence"}
    want_topk = bool(topk_consumers & enabled)
    topk_args = ["--top_k", str(args.top_k)] if want_topk else []

    for tag, fa in datasets:
        is_train = tag == "train"
        steps.append(Step(f"motif_{tag}", "motif_search.sh", ["--fasta_path", fa],
                          out_motif(fa), tool="motif"))
        steps.append(Step(f"motif_pair_{tag}", "motif_pair_distance.sh",
                          ["--fasta_path", fa], out_motif_pair(fa), tool="motif_pair"))
        steps.append(Step(f"esm_{tag}", "esm_embedding.sh", ["--fasta_path", fa],
                          out_esm(fa), tool="esm"))
        steps.append(Step(f"esm_ppl_{tag}", "esm_pseudo_perplexity.sh",
                          ["--fasta_path", fa], out_esm_ppl(fa), tool="esm_ppl"))
        steps.append(Step(f"maxid_self_{tag}", "max_sequence_identity.sh",
                          ["--fasta_path", fa] + (["--train"] if is_train else []),
                          out_maxid_self(fa), tool="maxid_self"))
        # Fast LOCAL (MMseqs2) sequence search within the dataset (self mode -> no
        # --train_path; the tool excludes each query's self-hit). Within-set novelty
        # counterpart of maxid_self. Explicit --save_path -> _self.csv so it does not
        # collide with the gen-vs-train step's default <fa>_local_sequence_search.csv.
        steps.append(Step(f"local_search_{tag}", "local_sequence_search.sh",
                          ["--fasta_path", fa, "--save_path", out_local_search_self(fa)],
                          out_local_search_self(fa), tool="local_sequence_search"))
        steps.append(Step(f"mindist_self_{tag}", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(fa)] + (["--train"] if is_train else []),
                          out_mindist_self(fa), tool="mindist_self", deps=[f"esm_{tag}"]))
        steps.append(Step(f"soluprot_{tag}", "soluprot.sh", ["--fasta_path", fa],
                          out_soluprot(fa), tool="soluprot"))
        steps.append(Step(f"ee_seq_{tag}", "enzyme_explorer_sequence_only.sh",
                          ["--fasta_path", fa], out_ee_seq(fa), tool="ee_seq"))

    # Broad "what-else-is-it-like" sequence search vs Swiss-Prot (annotated -> each hit
    # labeled TPS/non-TPS). Gen-only: it evaluates the generated designs (real train TPS
    # are trivially TPS hits).
    steps.append(Step("swissprot_search_gen", "swissprot_search.sh",
                      ["--fasta_path", gen], out_swissprot_search(gen), tool="swissprot_search"))

    if train:
        train_emb = args.train_embeddings_path or out_esm(train)
        # max_sequence_identity stays as the plain GLOBAL novelty metric (no --top_k):
        # the fast local_sequence_search now feeds the sequence-space top-k instead.
        steps.append(Step("maxid_gen_vs_train", "max_sequence_identity.sh",
                          ["--fasta_path", gen, "--train_path", train], out_maxid(gen),
                          tool="maxid_gen_vs_train"))
        # Fast LOCAL (MMseqs2) gen-vs-train search. Its _topk.csv (query_id,rank,
        # neighbour_id,score=identity%) is the SEQUENCE-space feeder for knn + sdr.
        steps.append(Step("local_search_gen_vs_train", "local_sequence_search.sh",
                          ["--fasta_path", gen, "--train_path", train] + topk_args,
                          out_local_search(gen), tool="local_sequence_search"))
        steps.append(Step("mindist_gen_vs_train", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(gen), "--train_embeddings_path", train_emb]
                          + topk_args,
                          out_mindist(gen), tool="mindist_gen_vs_train", deps=["esm_gen", "esm_train"]))

    # Structure branch: pLDDT (folding confidence) + foldseek structural identity to the
    # nearest known TPS. Both are keyed by the structures dir, run only if --structs_dir
    # is given. They have no SLURM deps (structures are produced out-of-band — by the
    # AlphaFold fan-out, not yet ported here) but plots picks them up via the snapshot.
    if structs:
        steps.append(Step("plddt_gen", "plddt.sh", ["--structs_dir", structs],
                          out_plddt(structs), tool="plddt"))
        steps.append(Step("motif_struct_gen", "motif_structural_distance.sh",
                          ["--structs_dir", structs], out_motif_struct(structs), tool="motif_struct"))
        # Active-site geometry: carboxylate-cage convergence (apo-robust). The
        # constellation-RMSD columns need reference templates present in the structs
        # dir (pass --templates <ids>); omitted here until a curated reference set is wired.
        steps.append(Step("active_site_geom_gen", "active_site_geometry.sh",
                          ["--structs_dir", structs], out_active_site_geom(structs),
                          tool="active_site_geom"))
        # Aromatic / cation-π pocket lining (carbocation-stabilization proxy).
        steps.append(Step("aromatic_lining_gen", "aromatic_lining.sh",
                          ["--structs_dir", structs], out_aromatic_lining(structs),
                          tool="aromatic_lining"))
        # Diphosphate-sensor basic residues (Arg/Lys + RY pair) at the metal site.
        steps.append(Step("diphosphate_sensor_gen", "diphosphate_sensor.sh",
                          ["--structs_dir", structs], out_diphosphate_sensor(structs),
                          tool="diphosphate_sensor"))
        # Radius of gyration / compactness: raw geometric shape numbers (Rg, asphericity,
        # principal radii) over the Cα atoms; no expected-Rg band (compared downstream).
        steps.append(Step("radius_of_gyration_gen", "radius_of_gyration.sh",
                          ["--structs_dir", structs], out_radius_of_gyration(structs),
                          tool="radius_of_gyration"))
        # Active-site pocket descriptors: fpocket geometry (catalytic-pocket volume,
        # hydrophobicity, enclosure) + P2Rank ML ligandability, anchored on the
        # carboxylate-cage metal point. Raw numbers; band from the reference-stats pipeline.
        steps.append(Step("pocket_descriptors_gen", "pocket_descriptors.sh",
                          ["--structs_dir", structs], out_pocket_descriptors(structs),
                          tool="pocket_descriptors"))
        # TPS structural-domain composition via EnzymeExplorer's CPU domain detector.
        steps.append(Step("domain_composition_gen", "domain_composition.sh",
                          ["--structs_dir", structs], out_domain_composition(structs),
                          tool="domain_composition"))
        # Aggrescan3D structure-based aggregation propensity (expressibility signal).
        steps.append(Step("aggregation_gen", "aggregation.sh",
                          ["--structs_dir", structs], out_aggregation(structs), tool="aggregation"))
        # Broad structural search vs AlphaFold-Swiss-Prot (annotated -> hit TPS/non-TPS).
        steps.append(Step("foldseek_swissprot_gen", "foldseek_swissprot_search.sh",
                          ["--structs_dir", structs], out_foldseek_swissprot(structs),
                          tool="foldseek_swissprot"))
        # ProteinMPNN sequence-likelihood (NLL) of the design's own sequence given its fold.
        steps.append(Step("proteinmpnn_gen", "proteinmpnn_score.sh",
                          ["--structs_dir", structs], out_proteinmpnn(structs), tool="proteinmpnn"))
        # Self-consistency scRMSD (ProteinMPNN -> ESMFold refold -> RMSD). HEAVY
        # (~1-2.5 min/structure x num_seqs GPU), so opt-in (config default off;
        # enable via --include self_consistency or the --self_consistency alias).
        steps.append(Step("self_consistency_gen", "self_consistency.sh",
                          ["--structs_dir", structs, "--num_seqs", str(args.sc_num_seqs)],
                          out_self_consistency(structs), tool="self_consistency"))
        if known_structs:
            steps.append(Step("structural_identity_gen", "structural_identity.sh",
                              ["--structs_dir", structs, "--known_structs_dir", known_structs]
                              + topk_args,
                              out_structural_identity(structs), tool="structural_identity"))
        else:
            print("[note] --structs_dir without --known_structs_dir: skipping "
                  "structural_identity (needs a known-TPS reference structure dir).")
        # PAE-consuming metrics need the per-structure PAE matrices saved at fold time
        # (<ID>_pae.npz, produced by ESMFold / the AF3 extract_pae step). Gated on
        # --pae_dir, like structural_identity is gated on --known_structs_dir.
        if pae_dir:
            steps.append(Step("global_confidence_gen", "global_confidence.sh",
                              ["--structs_dir", structs, "--pae_dir", pae_dir],
                              out_global_confidence(structs), tool="global_confidence"))
            steps.append(Step("interdomain_pae_gen", "interdomain_pae.sh",
                              ["--structs_dir", structs, "--pae_dir", pae_dir],
                              out_interdomain_pae(structs), tool="interdomain_pae"))
        else:
            print("[note] --structs_dir without --pae_dir: skipping global_confidence + "
                  "interdomain_pae (need <ID>_pae.npz from a PAE-saving fold).")

    # --- Dependency-heavy consumers: k-NN label transfer + SDR divergence --------- #
    # Both read the three feeders' *_topk.csv. They are gated on the feeders actually
    # existing (train -> sequence/embedding feeders; structs+known_structs -> structural
    # feeder) and skip cleanly with a [note] otherwise, like structural_identity does
    # without --known_structs_dir.
    #
    # k-NN: ensembled coarse-label transfer over the three similarity spaces. predict
    # mode needs ALL THREE feeders present, plus the label_file + calibration artifacts.
    # The committed calibration covers all three spaces (sequence via the fast MMseqs2
    # local_sequence_search). knn_calibrated_spaces() forwards only the calibrated spaces,
    # so a 2-space calibration would still work (sequence space abstained).
    if "knn" in enabled:
        knn_ok = bool(train and structs and known_structs)
        if not knn_ok:
            print("[note] knn: skipping label transfer (needs --train_path AND "
                  "--structs_dir AND --known_structs_dir for the three feeders).")
        elif not os.path.isfile(args.knn_label_file):
            print(f"[note] knn: skipping (label file not found: {args.knn_label_file}).")
        elif not os.path.isfile(args.knn_calibration):
            print(f"[note] knn: skipping (calibration not found: {args.knn_calibration}).")
        else:
            # Only forward the topk for spaces the calibration actually covers:
            # transfer_labels indexes calibration["spaces"][space] for every passed
            # space, so an uncalibrated one (e.g. sequence while only embedding+
            # structural are calibrated) would KeyError. space -> (flag, topk, feeder).
            cal_spaces = knn_calibrated_spaces(args.knn_calibration)
            space_specs = {
                "sequence":   ("--sequence_topk",   out_local_search_topk(gen),
                               "local_search_gen_vs_train"),
                "embedding":  ("--embedding_topk",  out_mindist_topk(out_esm(gen)),
                               "mindist_gen_vs_train"),
                "structural": ("--structural_topk", out_structural_identity_topk(structs),
                               "structural_identity_gen"),
            }
            knn_args = ["predict"]
            knn_deps: List[str] = []
            for space in ("sequence", "embedding", "structural"):
                if space not in cal_spaces:
                    continue
                flag, topk_path, feeder = space_specs[space]
                knn_args += [flag, topk_path]
                knn_deps.append(feeder)
            if not knn_deps:
                print(f"[note] knn: skipping (calibration covers no known space: "
                      f"{sorted(cal_spaces)}).")
            else:
                if "sequence" not in cal_spaces:
                    print("[note] knn: sequence space not yet calibrated -> using "
                          f"{sorted(s for s in cal_spaces if s in space_specs)} only.")
                knn_args += ["--label_file", args.knn_label_file,
                             "--calibration", args.knn_calibration,
                             "--out", out_knn(gen)]
                steps.append(Step("knn_gen", "knn_label_transfer.sh", knn_args,
                                  out_knn(gen), tool="knn", deps=knn_deps))

    # SDR divergence: needs the structural feeder (structs + known_structs); the sequence
    # feeder (--sequence_topk) is an optional fallback added only when train is present.
    if "sdr_divergence" in enabled:
        if not (structs and known_structs):
            print("[note] sdr_divergence: skipping (needs --structs_dir AND "
                  "--known_structs_dir).")
        else:
            sdr_args = ["--structs_dir", structs, "--known_structs_dir", known_structs,
                        "--structural_topk", out_structural_identity_topk(structs)]
            sdr_deps = ["structural_identity_gen"]
            if train:
                sdr_args += ["--sequence_topk", out_local_search_topk(gen)]
                sdr_deps.append("local_search_gen_vs_train")
            sdr_args += ["--save_path", out_sdr_divergence(structs)]
            steps.append(Step("sdr_divergence_gen", "sdr_divergence.sh", sdr_args,
                              out_sdr_divergence(structs), tool="sdr_divergence", deps=sdr_deps))

    # Filter to enabled tools BEFORE wiring plots, so plots only soft-depends on the
    # metric steps that actually survived (it never waits on a disabled tool).
    steps = [s for s in steps if s.tool in enabled]

    # plots: depend SOFTLY on every (enabled) metric (wait on whatever was submitted;
    # still run if some metrics were skipped/absent — the plots tool skips
    # missing-input targets). Only added if 'plots' itself is enabled.
    if "plots" in enabled:
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
                          tool="plots", soft_deps=metric_steps))
    return steps


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cluster", required=True, choices=sorted(CLUSTERS))
    p.add_argument("--fasta_path", required=True, help="Generated sequences FASTA.")
    p.add_argument("--train_path", default=None, help="Reference/train FASTA (enables comparisons).")
    p.add_argument("--structs_dir", default=None,
                   help="Generated structures dir (AF3 af_output or flat .pdb/.cif) -> "
                        "enables pLDDT and (with --known_structs_dir) structural-identity steps.")
    p.add_argument("--known_structs_dir", default=None,
                   help="Known-TPS reference structures dir for the foldseek structural-identity "
                        "metric (e.g. MARTS-DB train AFDB structures, or EE reference domains).")
    p.add_argument("--pae_dir", default=None,
                   help="Dir of per-structure PAE matrices (<ID>_pae.npz from a PAE-saving "
                        "ESMFold run or the AF3 extract_pae step) -> enables global_confidence "
                        "(pTM) and interdomain_pae.")
    p.add_argument("--train_structs_dir", default=None, help="(v2) train structures dir.")
    p.add_argument("--train_embeddings_path", default=None, help="Precomputed train embeddings CSV.")
    p.add_argument("--data_colors", nargs=2, default=["dodgerblue", "goldenrod"],
                   metavar=("TRAIN", "GEN"), help="Matplotlib colors for train/generated.")
    p.add_argument("--self_consistency", action="store_true",
                   help="Back-compat alias for '--include self_consistency': add the heavy "
                        "scRMSD self-consistency step (ProteinMPNN -> ESMFold refold -> RMSD; "
                        "~1-2.5 min/structure x num_seqs on GPU). Off by default.")
    p.add_argument("--sc_num_seqs", type=int, default=8,
                   help="ProteinMPNN sequences per structure for self_consistency (default 8).")
    # Single pipeline-wide top-k for the three feeders that emit *_topk.csv for the
    # knn + sdr_divergence consumers. Default comes from pipeline_tools.json
    # "_settings.top_k" (else the built-in DEFAULT_TOP_K); this overrides it.
    p.add_argument("--top_k", type=int, default=None,
                   help="Neighbours per query in the feeders' *_topk.csv (consumed by knn + "
                        f"sdr_divergence). Default: pipeline_tools.json _settings.top_k or "
                        f"{DEFAULT_TOP_K}. Only emitted when knn or sdr_divergence is enabled.")
    p.add_argument("--knn_label_file", default=DEFAULT_KNN_LABEL_FILE,
                   help="CSV mapping reference_id,label for k-NN label transfer "
                        f"(default {os.path.relpath(DEFAULT_KNN_LABEL_FILE, REPO)}).")
    p.add_argument("--knn_calibration", default=DEFAULT_KNN_CALIBRATION,
                   help="k-NN calibration JSON from `calibrate` "
                        f"(default {os.path.relpath(DEFAULT_KNN_CALIBRATION, REPO)}; "
                        "sequence+embedding+structural spaces).")
    # Config-driven tool selection (precedence: config defaults -> only -> include -> exclude).
    p.add_argument("--only", default=None, metavar="A,B,...",
                   help="Run ONLY these tool keys (comma-separated); plots is kept on unless "
                        "also excluded. See --list-tools for valid keys.")
    p.add_argument("--include", default=None, metavar="A,B,...",
                   help="Force-enable these tool keys on top of the config defaults.")
    p.add_argument("--exclude", default=None, metavar="A,B,...",
                   help="Force-disable these tool keys.")
    p.add_argument("--list-tools", action="store_true", dest="list_tools",
                   help="Print each tool key (default/branch/description) and exit.")
    p.add_argument("--dry-run", action="store_true", help="Print the submission plan; don't submit.")
    args = p.parse_args()

    catalog = load_tools_catalog()
    if args.list_tools:
        width = max(len(k) for k in catalog)
        print(f"Tools catalog ({TOOLS_CONFIG if os.path.isfile(TOOLS_CONFIG) else 'built-in defaults'}):")
        print(f"  {'KEY'.ljust(width)}  DEFAULT  BRANCH      DESCRIPTION")
        for key in catalog:
            v = catalog[key]
            on = "on " if v.get("default") else "off"
            print(f"  {key.ljust(width)}  {on}      {v.get('branch','?'):<10}  {v.get('description','')}")
        return

    args.fasta_path = os.path.abspath(args.fasta_path)
    for opt in ("train_path", "structs_dir", "known_structs_dir", "pae_dir", "train_structs_dir",
                "train_embeddings_path", "knn_label_file", "knn_calibration"):
        if getattr(args, opt):
            setattr(args, opt, os.path.abspath(getattr(args, opt)))

    # Resolve the single pipeline top-k: CLI --top_k wins, else JSON _settings/DEFAULT.
    if args.top_k is None:
        args.top_k = load_top_k()

    ts = datetime.now().strftime("%Y_%m_%d_%H%M%S")
    log_dir = os.path.join(REPO, "logs", f"run_eval_pipeline-{ts}")
    if not args.dry_run:
        os.makedirs(log_dir, exist_ok=True)

    enabled = resolve_enabled_tools(catalog, args)
    steps = build_steps(args, enabled)
    print(f"=== tps_eval pipeline on '{args.cluster}' — {len(steps)} steps; logs: {log_dir} ===")
    Engine(args.cluster, log_dir, dry_run=args.dry_run).run(steps)


if __name__ == "__main__":
    main()
